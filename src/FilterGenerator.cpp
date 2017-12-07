#ifndef VAR_FILTERGENERATOR_CC
#define VAR_FILTERGENERATOR_CC

#include "FilterGenerator.h"
#include "wavelet2d.h"
#include "AudioFile.h"
#include "gnuplot-iostream.h"

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <utility>


using namespace std;


FilterGenerator::FilterGenerator(vector<string> samples){
    level = 4;
    test("../data/in/positive.wav");
    vector<vector<double>> s_list;
    s_list = open_samples(samples);
    if (s_list.size() > 0){
        init_filter(s_list);
    }
    flag_stack.clear();
    length_stack.clear();
}


bool FilterGenerator::filter(string in_file, string out_file){
    vector<vector<double>> frame_list, frame_tree;
    vector<double> ref_frame, frame, d_e_list, d_e;
    AudioFile<double> af;
    if (af.load(in_file))
        init_frames(af.samples[0], frame_list);
    else
        return false;
    update_filter(af.samples[0]);
    for (int i = 0; i < (int)frame_list.size(); i++){
        cout << (double)i / (int)frame_list.size() * 100.0 << endl;
        frame = frame_list[i];
        ref_frame = frame;
        wpt_decompose(frame, frame_tree);
        apply_filter(frame_tree);
        iwpt_recompose(frame_tree, frame);
        d_e = energy_increase(ref_frame, frame);
        for (int j = 0; j < (int)d_e.size(); j++){
            d_e_list.push_back(d_e[j]);
        }
        frame_tree.clear();
        frame.clear();
        ref_frame.clear();
        d_e.clear();
    }
    plot_energy(d_e_list);
    thresshold_energy(d_e_list);
    suppress_noise(af.samples[0], d_e_list);
    cin.get();
    return af.save(out_file);
}


void FilterGenerator::update_filter(vector<double> n){
    vector<vector<double>> tree;
    wpt_decompose(n, tree);
    flag_stack.clear();
    length_stack.clear();
    e_noise = extract_energy(tree);
    H.resize(e_sig.size());
    for (int i = 0; i < (int)e_sig.size(); i++){
        H[i] = e_sig[i] / e_noise[i];
    }
}


void FilterGenerator::suppress_noise(vector<double>& n, vector<double> e){
    for (int i = 0; i < (int)e.size(); i++){
        if (e[i] == 1){
            n[i] /= 0.5;
        }
    }
}


void FilterGenerator::plot_energy(vector<double> e){
    Gnuplot gp;
    vector<pair<int, double>> plot_obj;
    for (int i = 0; i < (int)e.size(); i++){
        plot_obj.push_back(make_pair(i, e[i]));
    }
    gp << "plot '-'\n";
    gp.send1d(plot_obj);
}


void FilterGenerator::thresshold_energy(vector<double>& e){
    for (int i = 0; i < (int)e.size(); i++){
        if (e[i] < 0.045)
            e[i] = 0;
        else
            e[i] = 1;
    }
}


void FilterGenerator::apply_filter(vector<vector<double>>& frame_tree){
    for (int i = 0; i < (int)frame_tree.size(); i++){
        for (int j = 0; j < (int)frame_tree[i].size(); j++){
            frame_tree[i][j] *= H[i];
        }
    }
}


vector<double> FilterGenerator::normalize_energy(vector<double> sig){
    vector<double> e;
    e.resize(frame_size);
    double t_e = 0.0;
    for (int i = 0; i < (int)sig.size(); i++){
        t_e += abs(sig[i] * sig[i]);
    }
    for (int i = 0; i < (int)sig.size(); i++){
        e[i] = abs(sig[i] * sig[i]) / t_e;
    }
    return e;
}


vector<double> FilterGenerator::energy_increase(vector<double> x, vector<double> y){
    vector<double> e_in, e_out, d_e;
    double t;
    e_in = normalize_energy(x);
    e_out = normalize_energy(y);
    d_e.resize(frame_size);
    for (int i = 0; i < (int)d_e.size(); i++){
        t = 0.0;
        for (int j = i; j < frame_size; j++){
            t += e_out[j] - e_in[j];
        }
        d_e[i] = t;
    }
    return d_e;
}


void FilterGenerator::init_frames(vector<double> sig, vector<vector<double>>& frame_list){
    vector<double> frame;
    int i = 0;
    siglen = (int)sig.size();
    skip = frame_size;
    while (i < siglen){
        if (i + frame_size < siglen){ 
            for (int j = 0; j < frame_size; j++){
                frame.push_back(sig[i + j]);
            }
            frame_list.push_back(frame);
            frame.clear();
        }
        i += skip;
    }
}


//open files and load into array
vector<vector<double>> FilterGenerator::open_samples(vector<string> file_name){
    vector<vector<double>> s_list;
    AudioFile<double> af;
    for (int i = 0; i < (int)file_name.size(); i++){
        if (af.load(file_name[i])){
            s_list.push_back(af.samples[0]);
        } else {
            cout << file_name[i] << " invalid filename!" << endl;
        }
    }
    //TODO: add padding to make all samples same size
    frame_size = (int)s_list[0].size();
    return s_list;
}


void FilterGenerator::init_filter(vector<vector<double>> s_list){
    vector<vector<vector<double>>> wpt_tree_list;
    vector<vector<double>> wpt_avg_tree;
    wpt_tree_list.resize(s_list.size());
    //decompose all signals
    for (int i = 0; i < (int)s_list.size(); i++){
        wpt_decompose(s_list[i], wpt_tree_list[i]);
    }
    //calculate average wpt
    avg_tree_list(wpt_tree_list, wpt_avg_tree);

    e_sig = extract_energy(wpt_avg_tree);
}


vector<double> FilterGenerator::extract_energy(vector<vector<double>> wpt_avg_tree){
    vector<double> e;
    e.resize(wpt_avg_tree.size());
    double t;
    for (int i = 0; i < (int)wpt_avg_tree.size(); i++){
        e[i] = 0;
        for (int j = 0; j < (int)wpt_avg_tree[i].size(); j++){
            t = abs(wpt_avg_tree[i][j]);
            e[i] += t * t;
        }
        e[i] *= 1.0 / (int)wpt_avg_tree[i].size();
    }
    return e;
}


void FilterGenerator::avg_tree_list(vector<vector<vector<double>>> wpt_tree_list, 
        vector<vector<double>>& wpt_avg_tree){
    wpt_avg_tree = wpt_tree_list[2];
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
        orig.resize(orig.size()/8);
        wpt = orig;

        wpt_decompose(wpt, wpt_tree);
        iwpt_recompose(wpt_tree, wpt);
        

        //test
        bool retval = true;



        cout << wpt.size() << ", " << orig.size() << endl;
        double diff = 0;
        double origtotal = 0;
        double newtotal = 0;
            cout << "same size" << endl;
        for (int i = 0; i < (int)wpt.size(); i++){
            if (wpt[i] != orig[i]){
                retval = false;
                diff += abs(wpt[i] - orig[i]);
            }
            origtotal += abs(orig[i]);
            newtotal += abs(wpt[i]);
        }
        cout << "diff: " << diff << endl;
        cout << "orig total: " << origtotal << endl;
        cout << "new total: " << newtotal << endl;
        return retval;
    }
    return false;
}


void FilterGenerator::wpt_decompose(vector<double> in, vector<vector<double>>& out){
    vector<vector<vector<double>>> cwt;
    vector<double> temp, flag;
    vector<int> length;

    cwt.resize(level);
    for (int i = 0; i < level; i++){
        cwt[i].resize(1<<i);
    }

    //root node tree
    cwt[0][0] = in;

    cout << in.size() << endl;
    for (int i = 0; i < (int)cwt.size()-1; i++){
        for (int j = 0; j < (int)cwt[i].size(); j++){
            dwt(cwt[i][j], 1, "db3", temp, flag, length);
            flag_stack.push_back(flag);
            length_stack.push_back(length);
            for (int k = 0; k < (int)temp.size(); k++){
                if (length[0] > k){
                    cwt[i+1][j*2].push_back(temp[k]);
                } else {
                    cwt[i+1][j*2+1].push_back(temp[k]);  
                }
            }
            temp.clear();
            flag.clear();
            length.clear();
        }
    }
    out = cwt[level - 1];
}


void FilterGenerator::iwpt_recompose(vector<vector<double>> in, vector<double>& out){
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
            for (int k = 0; k < (int)cwt[i+1][j*2].size(); k++){
                temp.push_back(cwt[i+1][j*2][k]);
            }
            for (int k = 0; k < (int)cwt[i+1][j*2+1].size(); k++){
                temp.push_back(cwt[i+1][j*2+1][k]);
            }
            idwt(temp, flag_stack[s], "db3", cwt[i][j], length_stack[s]);
            temp.clear();
            s--;
        }
    }
    out = cwt[0][0];
    flag_stack.clear();
    length_stack.clear();
}



#endif
