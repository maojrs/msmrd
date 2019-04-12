//
// Created by maojrs on 2/6/19.
//

#include "markovModels/msmrdMarkovModel.hpp"

namespace msmrd{
    using ctmsm = msmrd::continuousTimeMarkovStateModel;
    /**
    * Implementation of Markov model class specializes for MSM/RD.
    */
    msmrdMarkovModel::msmrdMarkovModel(unsigned int numBoundStates, unsigned int numTransitionStates,
                                       long seed, std::map<std::string, float> &rateDictionary)
            : numBoundStates(numBoundStates), numTransitionStates(numTransitionStates),
              seed(seed), rateDictionary(rateDictionary) {
        randg.setSeed(seed);
    };

    /* Computes transition time and end state starting from a given transition step. The end states correspond
     * to one of the bound states (indexing of bound states begins in 1)*/
    std::tuple<double, int> msmrdMarkovModel::computeTransition(int transitionState) {
        std::vector<float> rates(numBoundStates);
        double transitionTime;
        int endState = 0;

        // Get rates to all bound states from dictionary
        std::string key;
        for (int i = 0; i < numBoundStates; i++) {
            key = std::to_string(transitionState) + "->b" + std::to_string(i+1);
            rates[i] = getRate(key);
        }

        // Calculate lambda0  and ratescumsum for SSA/Gillespie algorithm
        double lambda0 = 0;
        for(const auto& rate: rates) {
            lambda0 += rate;
        }

        // Calculate time for transition
        double r1 = randg.uniformRange(0, 1);
        transitionTime = std::log(1.0/r1)/lambda0;

        // Calculate end state
        double r2lam0 = randg.uniformRange(0, 1)*lambda0;
        double cumsum = 0;
        for (int i=0; i< numBoundStates; i++){
            cumsum += rates[i];
            if (r2lam0 <= cumsum) {
                endState = i+1;
                break;
            }
        }

        return std::make_tuple(transitionTime, endState);
    }


    /* Getter functions for dictionary, a bit too long to leave in headers.
     * Better than simply calling dicionary[key] since that might create a new key with empty
     * value when called.*/

    /* @param key is in the form of "state1->state2", e.g. "b1->14" or "32->b2",
     * where the states starting with "b" are bound states and the rest correspond to transition states*/
    float msmrdMarkovModel::getRate(std::string key) {
        auto search = rateDictionary.find(key);
        if (search != rateDictionary.end()) {
            return search->second;
        } else {
            throw std::range_error("Key does not exist in dictionary");
        }
    }

}
