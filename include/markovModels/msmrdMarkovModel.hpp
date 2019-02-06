//
// Created by maojrs on 2/6/19.
//

#pragma once
#include <map>
#include "markovModels/continuousTimeMarkovModel.hpp"

namespace msmrd {
    using ctmsm = continuousTimeMarkovStateModel;
    /*
     * Markov state model class used for MSM/RD algorithm. It is based on the continuous-time
     * Markov state model class. However, it has additional functionality that is required for the
     * coupling.
     */
    class msmrdMarkovModel : public ctmsm {
    protected:


        std::map<std::string, float> rateDictionary;
    public:
        // Inherit constructor from parent and add an additional constructor to read rates.
        //using ctmsm::ctmsm;

        msmrdMarkovModel(int msmid, std::vector<std::vector<double>> &tempmatrix, long seed,
                         std::map<std::string, float> &rateDictionary);
    };

}


