//
// Created by maojrs on 2/6/19.
//

#include "markovModels/msmrdMarkovModel.hpp"

namespace msmrd{
    /**
     * Implementation of Markov model class specializes for MSM/RD.
     */
    msmrdMarkovModel::msmrdMarkovModel(unsigned int numBoundStates,
                                                 unsigned int numTransitionStates, long seed,
                                                 std::map<std::string, float> &rateDictionary)
            : numBoundStates(numBoundStates), numTransitionStates(numTransitionStates), seed(seed),
              rateDictionary(rateDictionary) {
        randg.setSeed(seed);
        Dboundlist.resize(numBoundStates);
        DboundRotlist.resize(numBoundStates);
    };

    /* Computes transition time and end state starting from a given transition step. The end states correspond
     * to one of the bound states (indexing of bound states begins in 1)*/
    std::tuple<double, int> msmrdMarkovModel::computeTransition2BoundState(int transitionState) {
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
        for(const auto &rate: rates) {
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


    /* Computes transition time and end state starting from a given bound state. The end states correspond
     * to either another bound state or a transition state (indexing of bound states begins in 1)*/
    std::tuple<double, int> msmrdMarkovModel::computeTransitionFromBoundState(int boundState) {
        std::vector<float> rates(maxNumberBoundStates + numTransitionStates);
        double transitionTime;
        int index0 = maxNumberBoundStates;
        int endState = 0;

        // Get rates to all bound states from dictionary
        std::string key;
        for (int i = 0; i < maxNumberBoundStates; i++) {
            key = "b" + std::to_string(boundState) + "->b" + std::to_string(i + 1);
            rates[i] = getRate(key); // non-existent keys yield rate zero
        }
        // Get rates to all transition states from dictionary
        for (int i = 0; i < numTransitionStates; i++) {
            key = "b" + std::to_string(boundState) + "->" + std::to_string(i + 1); //index0 not required here
            rates[index0 + i] = getRate(key);
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
        for (int i=0; i < maxNumberBoundStates + numTransitionStates; i++){
            cumsum += rates[i];
            if (r2lam0 <= cumsum) {
                endState = i + 1;
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
            return 0.0;
            //throw std::range_error("Key does not exist in dictionary");
        }
    }

    /* These functions set the diffusion coefficients for bound and unbound states (A and B).
     * Note size of vectors must match numBoundStates. This functions NEEDS
     * to be called to fill in diffusion coefficients. */

    void msmrdMarkovModel::setDbound(std::vector<double> &D, std::vector<double> &Drot) {
        if ( (D.size() != numBoundStates ) or (Drot.size() != numBoundStates) ) {
            std::__throw_range_error("Vectors of diffusion coefficients must match number of bound states");
        }
        Dboundlist = D;
        DboundRotlist = Drot;
    }

}
