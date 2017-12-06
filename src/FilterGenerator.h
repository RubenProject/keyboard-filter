#ifndef VAR_FILTERGENERATOR_H
#define VAR_FILTERGENERATOR_H

#include "wavelet2d.h"
#include "AudioFile.h"

#include <string>
#include <vector>

using namespace std;


class FilterGenerator {
    public:
        FilterGenerator(vector<string> samples);
        bool test(string file_name);
        bool filter(string in_file, string out_file);
    private:
        void wpt_decompose(vector<double> in, vector<vector<double>>& out);
        void iwpt_recompose(vector<vector<double>> in, vector<double>& out);
        int make_p2(vector<double>& samples);

        vector<vector<double>> open_samples(vector<string> file_name);
        void init_filter(vector<vector<double>> s_list);
        void avg_tree_list(vector<vector<vector<double>>> wpt_tree_list, vector<vector<double>>& wpt_avg_tree);
        vector<double> extract_energy(vector<vector<double>> wpt_avg_tree);

        void init_frames(vector<double> s, vector<vector<double>>& frame_list);
        void apply_filter(vector<vector<double>>& frame_tree);
        double energy_increase(vector<double> ref, vector<double> f);

        int level;
        vector<vector<double>> flag_stack;
        vector<vector<int>> length_stack;

        vector<double> H;
        int siglen;
        int frame_size;
        int skip;
};


#endif
