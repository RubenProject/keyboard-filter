#include "FilterGenerator.h"

#include <iostream>
#include <string>
#include <vector>




int main(){
    std::vector<std::string> sample_name = {"../data/samples/1.wav",
                                            "../data/samples/2.wav",
                                            "../data/samples/3.wav",
                                            "../data/samples/4.wav",
                                            "../data/samples/5.wav",
                                            "../data/samples/6.wav"};

    FilterGenerator FG(sample_name);
    //FG.filter("../data/in/labeled_negative.wav", "../data/out/labeled_negative.wav");
    FG.filter("../data/in/labeled_positive.wav", "../data/out/labeled_positive.wav");

    return 0;
}
