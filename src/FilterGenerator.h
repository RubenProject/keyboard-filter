#ifndef VAR_FILTERGENERATOR_H
#define VAR_FILTERGENERATOR_H

#include "../wavelib/header/wauxlib.h"
#include "../wavelib/header/wavelib.h"
#include "../AudioFile/AudioFile.h"

#include <string>
#include <vector>

class FilterGenerator {
    public:
        FilterGenerator(std::vector<std::string> samples);
        ~FilterGenerator();
    private:
        void extract_energy(AudioFile<double> af);

};


#endif
