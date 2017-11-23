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



FilterGenerator::FilterGenerator(std::vector<std::string> samples){
    AudioFile<double> af;
    for (int i = 0; i < (int)samples.size(); i++){
        if (af.load(samples[i])){
            extract_energy(af);
        }
    }
}

FilterGenerator::~FilterGenerator(){
}


void FilterGenerator::extract_energy(AudioFile<double> af){
    std::cout << "attempting transform" << std::endl;
    wave_object obj;
    wt_object wt;
    double *in, *out;
    int N, i, J;
    
    char* name = "haar";
    obj = wave_init(name);

    N = (int)af.samples[0].size();

    in = &af.samples[0][0];

    J = 2; //2th level decomposition

    wt = wt_init(obj, "dwt", N, J);
    dwt(wt, in);
    wt_summary(wt);

}

#endif
