#include "FilterGenerator.h"

#include <iostream>
#include <string>
#include <vector>




int main(){
    std::vector<std::string> sample_name = {"../data/samples/key_sharp.wav",
                                            "../data/samples/key1.wav",
                                            "../data/samples/key2.wav"};

    FilterGenerator FG(sample_name);
    FG.filter("../data/in/p_real.wav", "../data/out/p_real.wav");
    //FG.filter("../data/in/positive.wav", "../data/out/positive.wav");
    //FG.filter("../data/in/negative.wav", "../data/out/negative.wav");

    



    return 0;
}
