#ifndef VAR_FILTERGENERATOR_CC
#define VAR_FILTERGENERATOR_CC

#include "FilterGenerator.h"

#include "../wavelib/header/wavelib.h"
#include "../wavelib/header/wauxlib.h"
#include "../AudioFile/AudioFile.h"

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <string>


using namespace std;



FilterGenerator::FilterGenerator(vector<string> samples){
    test(samples[0]);
    cout << "test successful" << endl;
    exit(0);
    vector<AudioFile<double>> af_list;
    level = 4;
    int sample_size;
    af_list = open_samples(samples, sample_size);

    vector<vector<double>> wpt_tree_set;
    wpt_tree_set.resize(af_list.size());

    cout << "Filter length: " << sample_size << endl;

    double* wpt_tree;
    for (int i = 0; i < (int)af_list.size(); i++){
        wpt_tree_set[i].resize(sample_size);
        wpt_tree = (double*)malloc(sample_size * sizeof(double));
        for (int j = 0; j < sample_size; j++){
            wpt_tree[j] = af_list[i].samples[0][j];
        }
        wpt_decompose(wpt_tree, sample_size , level);
        for (int j = 0; j < sample_size; j++){
            wpt_tree_set[i][j] = wpt_tree[j];
        }
        free(wpt_tree);
    }
    
    vector<double> wpt_tree_avg;
    wpt_tree_avg = avg_tree(wpt_tree_set, sample_size);
    //there exist 2^level subbands
    int N = 1 << level;
    int M = sample_size / N;
    H = extract_energy(wpt_tree_avg, M, N);
    window_size = sample_size;
}


bool FilterGenerator::test(string file_name){
    AudioFile<double> test;
    int sample_size;
    int level;
    bool retval = true;
    if (test.load(file_name)){
        sample_size = test.samples[0].size();
        level = 4;
        double* wpt_tree = (double*)malloc((sample_size + 100) * sizeof(double));
        for (int i = 0; i < sample_size; i++){
            wpt_tree[i] = test.samples[0][i];
        }
        wpt_decompose(wpt_tree, sample_size, level);
        iwpt_recompose(wpt_tree, sample_size, level);
        for (int i = 0; i < sample_size; i++){
            if (wpt_tree[i] != test.samples[0][i]){
                cout << wpt_tree[i] << "!=" << test.samples[0][i] << endl;
                retval = false;
            }
        }
    }
    return retval;
}


void FilterGenerator::filter(string file_name){
    AudioFile<double> noise;
    double* frame = (double*)malloc(window_size * sizeof(double));
    if (noise.load(file_name)){
        int noise_length = (int)noise.samples[0].size();
        if (noise_length < 10 * window_size){
            cout << "Noise length too short" << endl;
            return;
        }
        for (int i = 0; i < noise_length; i += window_size){
            for (int j = 0; j < window_size; j++)
                frame[j] = noise.samples[0][i + j];
            wpt_decompose(frame, window_size, level);
            apply_filter(frame, window_size, 1 << level);
            //inverse WPT
            iwpt_recompose(frame, window_size, level);

            //meassure energy increase
            //normalize signals

            //apply thressholding
            //reduce the noise in signal
            //???
            //profit
        }


    } else {
        cout << "no such file" << endl;
    }
}


void FilterGenerator::apply_filter(double* frame, int M, int N){
    for (int i = 0; i < N; i++){
        for (int k = 0; k < M; k++){
            frame[i * M + k] *= H[i];
        }
    }
}


vector<double> FilterGenerator::extract_energy(vector<double> wpt_tree, int M, int N){
    vector<double> e;
    e.resize(N);
    double t;
    for (int i = 0; i < N; i++){
        e[i] = 0;
        for (int k = 0; k < M; k++){
            t = abs(wpt_tree[i * M + k]);
            e[i] += t * t;
        }
        e[i] *= 1 / M;

    }
    return e;
}


vector<AudioFile<double>> FilterGenerator::open_samples(vector<string>samples, int& sample_size){
    vector<AudioFile<double>> af_list;
    AudioFile<double> af;
    sample_size = 0;
    for (int i = 0; i < (int)samples.size(); i++){
        if (af.load(samples[i])){
            af_list.push_back(af);
            if (sample_size == 0 || (int)af.samples[0].size() < sample_size){
                sample_size = (int)af.samples[0].size();
            }
        }
    }
    return af_list;
} 

vector<double> FilterGenerator::avg_tree(vector<vector<double>> wpt_tree_set, int sample_size){
    vector<double> wpt_tree_avg;
    wpt_tree_avg.resize(sample_size);
    for (int i = 0; i < (int)wpt_tree_set.size(); i++){
        wpt_tree_avg[i] = 0;
        for (int j = 0; j < sample_size; j++){
            wpt_tree_avg[i] += wpt_tree_set[i][j];
        }
        wpt_tree_avg[i] /= sample_size;
    }
    return wpt_tree_avg;
}


void FilterGenerator::iwpt_recompose(double* wpt_tree, int sample_size, int level){
    wave_object obj;
    wt_object wt;

    double* out;
    out = (double*)malloc(100000 * sizeof(double));
    cout << "test" << endl;
    
    char* name = "db3";
    obj = wave_init(name);

    wt = wt_init(obj, "dwt", 100, 1);
    setDWTExtension(wt, "sym");
    wt->lenlength = 2;
    wt->length[0] = 160;
    wt->length[1] = 160;
    for (int i = 0; i < sample_size * 2; i++){
        wt->output[i] = wpt_tree[i];
        cout << wpt_tree[i] << endl;
    }

    
    wt_summary(wt);
    idwt(wt, out);
    
    wt_summary(wt);

    return; 




    for (int j = 0; j < 1 << (level - 1); j++){ 
        for (int i = 0; i < sample_size * 2; i++){
            wt->output[i] = wpt_tree[j * sample_size + i];
        }
        wt->siglength = sample_size;
        idwt(wt, out);
        wt_summary(wt);
        for (int i = 0; i < sample_size * 2; i++){
            wpt_tree[j * sample_size + i] = out[i];
        }
    }

}


void FilterGenerator::wpt_decompose(double* wpt_tree, int& sample_size, int level){
    wave_object obj;
    wtree_object wt;
    char *name = "db3";
    obj = wave_init(name);

    wt = wtree_init(obj, sample_size, level);
    setWTREEExtension(wt, "sym");

    wtree(wt, wpt_tree);
    wtree_summary(wt);

    int i = 0, node = 0;
    int length = getWTREENodelength(wt, level);

    while (i < sample_size){
        getWTREECoeffs(wt, level, node, &wpt_tree[i], length);
        node++;
        i += length;
    }

    sample_size = length;

    wave_free(obj);
    wtree_free(wt);
    
    
    

    //wave_object obj;
    //wt_object wt;
    
    //char* name = "haar";
    //obj = wave_init(name);
    //wt = wt_init(obj, "dwt", sample_size, 1);
    //setWTConv(wt, "fft");


    /*for (int i = 0; i < level; i++){
        for (int j = 0; j < 1 << i; j++){
            dwt(wt, &wpt_tree[sample_size * j]);
            //wt_summary(wt);
            setWTConv(wt, "direct");
            //copy in the new parts
            for (int k = 0; k < sample_size; k++){
                wpt_tree[sample_size * j + k] = wt->output[k];
            }
        }
        sample_size /= 2;
        wt->siglength /= 2;
    }
    */
}

#endif
