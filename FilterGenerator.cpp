#ifndef VAR_FILTERGENERATOR_CC
#define VAR_FILTERGENERATOR_CC

#include "FilterGenerator.h"
#include "wavelet2d.h"
#include "AudioFile.h"

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <string>


using namespace std;



FilterGenerator::FilterGenerator(vector<string> samples){
    level = 4;
    test(samples[0]);
}


int FilterGenerator::make_p2(vector<double>& samples){
    int n = (int)samples.size();
    while ((n & (n - 1)) != 0){
        n++;
    }
    samples.resize(n, 0.0);
    return n;
}


//test if wpt is close to almost lossless
bool FilterGenerator::test(string file_name){
    AudioFile<double> test;
    vector<vector<double>> wpt_tree;
    vector<double> wpt, orig;
    if (test.load(file_name)){
        orig = test.samples[0];
        cout << orig.size() << endl;
        make_p2(orig);
        cout << orig.size() << endl;
        wpt = orig;

        wpt_decompose(wpt, wpt_tree);
        iwpt_recompose(wpt_tree, wpt);



        //test
        cout << "----------------------"<< endl;
        bool retval = true;



        cout << wpt.size() << ", " << orig.size() << endl;
        double diff = 0;
        double origtotal = 0;
        double newtotal = 0;
        if (wpt.size() == orig.size())
            cout << "same size" << endl;
        for (int i = 0; i < (int)wpt.size(); i++){
            if (wpt[i] != orig[i]){
                retval = false;
                diff += abs(wpt[i] - orig[i]);
                origtotal += abs(orig[i]);
                newtotal += abs(wpt[i]);
            }
        }
        cout << "diff: " << diff << endl;
        cout << "orig total: " << origtotal << endl;
        cout << "new total: " << newtotal << endl;
        return retval;
    }
    return false;
}


void FilterGenerator::wpt_decompose(vector<double> in, vector<vector<double>>& out){
    cout << "-----------------------"<< endl;
    vector<vector<vector<double>>> cwt;
    vector<double> temp, flag;
    vector<int> length;

    cwt.resize(level);
    for (int i = 0; i < level; i++){
        cwt[i].resize(1<<i);
    }

    //root node tree
    cwt[0][0] = in;

    for (int i = 0; i < (int)cwt.size()-1; i++){
        for (int j = 0; j < (int)cwt[i].size(); j++){
            dwt_sym(cwt[i][j], 1, "db3", temp, flag, length);
            flag_stack.push_back(flag);
            length_stack.push_back(length);
            for (int k = 0; k < (int)temp.size(); k++){
                if (length[0] < k){
                    cwt[i+1][j*2].push_back(temp[i]);
                } else {
                    cwt[i+1][j*2+1].push_back(temp[i]);  
                }
            }
            cout << length[0] << endl;
            temp.clear();
            flag.clear();
            length.clear();
        }
    }
    out = cwt[level - 1];
}


void FilterGenerator::iwpt_recompose(vector<vector<double>> in, vector<double>& out){
    cout << "-----------------------"<< endl;
    vector<vector<vector<double>>> cwt;
    vector<double> temp, flag;
    vector<int> length;

    cwt.resize(level);
    for (int i = 0; i < level; i++){
        cwt[i].resize(1<<i);
    }

    //init leafs
    for (int i = 0; i < 1<<(level-1); i++){
        cwt[level - 1][i] = in[i];
    }

    int s = length_stack.size() - 1;

    //iterate tree backwards
    for (int i = (int)cwt.size() - 2; i >= 0; i--){
        for (int j = (int)cwt[i].size() - 1; j >= 0; j--){
            for (int k = 0; k < (int)(cwt[i+1][j*2].size() + cwt[i+1][j*2+1].size()); k++){
                if (k < (int)cwt[i+1][j*2].size()){
                    temp.push_back(cwt[i+1][j*2][k]);
                } else {
                    temp.push_back(cwt[i+1][j*2+1][k]);
                }
            }
            idwt_sym(temp, flag_stack[s], "db3", cwt[i][j], length_stack[s]);
            temp.clear();
            s--;
        }
    }
    out = cwt[0][0];
}



#endif
