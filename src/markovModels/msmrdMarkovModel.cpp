//
// Created by maojrs on 2/6/19.
//

#include "markovModels/msmrdMarkovModel.hpp"

namespace msmrd{
    /**
    * Implementation of Markov model class specializes for MSM/RD.
    * @param rateDict
    */
    msmrdMarkovModel::msmrdMarkovModel(int msmid, std::vector<std::vector<double>> &tempmatrix, long seed,
                                       std::map<std::string, float> &rateDictionary)
            : ctmsm::continuousTimeMarkovStateModel(msmid, tempmatrix, seed), rateDictionary(rateDictionary) {};

    // Need smart and fast way to load transition rates into class from dictionary.

}
