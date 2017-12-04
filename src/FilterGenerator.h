#ifndef VAR_FILTERGENERATOR_H
#define VAR_FILTERGENERATOR_H

#include "../wavelib/header/wauxlib.h"
#include "../wavelib/header/wavelib.h"
#include "../AudioFile/AudioFile.h"

#include <string>
#include <vector>

using namespace std;


class FilterGenerator {
    public:
        FilterGenerator(vector<string> samples);
        void filter(string file_name);
    private:
        void wpt_decompose(double* wpt_tree, int M, int N);
        vector<AudioFile<double>> open_samples(vector<string> samples, int& sample_size);
        vector<double> avg_tree(vector<vector<double>> wpt_tree_set, int sample_size);
        vector<double> extract_energy(vector<double> wpt_tree, int sample_size, int level);
        void apply_filter(double* frame, int M, int N);

        vector<double> H;
        int level;
        int window_size;
};


#endif
