#ifndef VAR_FILTERGENERATOR_H
#define VAR_FILTERGENERATOR_H

#include "wavelet2d.h"
#include "AudioFile.h"

#include <string>
#include <vector>

#define W_NAME "db3"
#define vec3d vector<vector<vector<double> > >
#define vec2d vector<vector<double> >
#define vec1d vector<double>
#define vec2i vector<vector<int> >
#define vec1i vector<int>


using namespace std;


class FilterGenerator {
    public:
        FilterGenerator(vector<string> samples);
        bool test(string file_name);
        bool filter(string in_file, string out_file);
        void print_stats();
    private:
        void wpt_decompose(vec1d in, vec2d& out);
        void iwpt_recompose(vec2d in, vec1d& out);

        vec1i calculate_labels(AudioFile<double> af);
        void calculate_error(vec1i l, vec1i d);

        vec2d open_samples(vector<string> file_name);
        void init_filter(vec2d s_list);
        void avg_energy(vec2d e);
        vec1d extract_energy(vec2d wpt_tree);

        void init_frames(vec1d s, vec2d& frame_list);
        void update_filter(vec1d n);
        void apply_filter(vec2d& frame_tree);
        vec1d energy_increase(vec1d ref, vec1d f);
        vec1d normalize_energy(vec1d sig);
        vec1i thresshold_energy(vec1d e, double t);
        void suppress_noise(vec1d& n, vec1d e);
        void plot_energy(vec1d e);


        int level;
        vec2d flag_stack;
        vec2i length_stack;

        vec1d e_sig;
        vec1d e_noise;
        vec1d H;
        int siglen;
        int frame_size;
        int skip;
};


#endif
