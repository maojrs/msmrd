//
// Created by maojrs on 2/6/19.
//

#include "markovModels/msmrdMarkovModel.hpp"

namespace msmrd{
    /**
    * Implementation of Markov model class specializes for MSM/RD.
    * @param rateDict
    */
    msmrdMarkovModel::msmrdMarkovModel(unsigned int nstates, unsigned int numDiscreteOrientations,
                                       long seed, std::map<std::string, float> &rateDictionary)
            : nstates(nstates), numDiscreteOrientations(numDiscreteOrientations),
              seed(seed), rateDictionary(rateDictionary) {
        randg.setSeed(seed);
        generateMarkovModels();
    };

    /* Generates several continuous times Markov models each representing a discretization of
     * the relative orientation */
    void msmrdMarkovModel::generateMarkovModels() {
        int numBoundStates = nstates - 1;
        int numMSMs = numBoundStates*(numBoundStates - 1) + (numDiscreteOrientations)*(numDiscreteOrientations + 1);
        ctmsmDictionary;

    }
}
