//
// Created by maojrs on 2/6/19.
//

#include "markovModels/msmrdMarkovModel.hpp"

namespace msmrd{
    using ctmsm = msmrd::continuousTimeMarkovStateModel;
    /**
    * Implementation of Markov model class specializes for MSM/RD.
    * @param nstates number of states in the msm (including unbound state)
    * @param numDiscreteOrientations number of Galerkin discrete orientations of one particle.
    * @param seed variable for random number generation (Note values of seed <= -1 correspond to random device)
    * @param rateDictionary dictionary relating transitions to its corresponding rates. The keys
    * are strings of the form "state1->state2". The bound states hava b before the number, and the transition
    * states are denoted with duple integers.
    */
    msmrdMarkovModel::msmrdMarkovModel(unsigned int nstates, unsigned int numDiscreteOrientations,
                                       long seed, std::map<std::string, float> &rateDictionary)
            : nstates(nstates), numDiscreteOrientations(numDiscreteOrientations),
              seed(seed), rateDictionary(rateDictionary) {
        generateMarkovModels();
    };

    /* Generates several continuous-time Markov models each representing a discretization of
     * the relative orientation */
    void msmrdMarkovModel::generateMarkovModels() {
        int numBoundStates = nstates - 1;
        // Create all MSMs corresponding to different orientations
        float rate;
        std::vector<std::vector<double>> transitionMatrix;
        int msmid = 0;
        for (int i = 1; i <= numDiscreteOrientations; i++){
            for (int j = i; j <= numDiscreteOrientations; j++){
                // Reset values of transition matrix to zero.
                transitionMatrix = std::vector<std::vector<double>> (nstates, std::vector<double> (nstates, 0));
                // Extract on rates from off-transition state ij to all bound states
                for (int k = 1; k <= numBoundStates; k++) {
                    auto key = std::to_string(i) + std::to_string(j) + "->b" + std::to_string(k);
                    auto search = rateDictionary.find(key);
                    if (search != rateDictionary.end()) {
                        rate = search->second;
                    } else {
                        rate = 0;
                    }
                    transitionMatrix[0][k] = rate;
                    transitionMatrix[0][0] -= rate;
                }
                // Create ctmsm and insert it into dictionary with key 'ij'.
                ctmsm ctmsm_ij = ctmsm(msmid, transitionMatrix, seed*(msmid+1));
                msmid += 1;
                auto newkey = std::to_string(i) + std::to_string(j);
                ctmsmDictionary.insert ( std::pair<std::string, ctmsm>(newkey, ctmsm_ij) );
            }
        }


//        // Extract off rates to unbound transition states
//        auto key2 = "b" + std::to_string(k) + "->" + std::to_string(i) + std::to_string(j);
//        auto search2 = rateDictionary.find(key2);
//        if (search2 != rateDictionary.end()) {
//            rate = search2->second;
//        } else {
//            rate = 0;
//        }
//        transitionMatrix[k][0] = rate;
//        transitionMatrix[k][0] -= rate;

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

    /* @param key is in the form of "transition_region", e.g. "14" or "32", which correspond to the two closest
     * orientational discrete states of each of the two particles (touched by a line)
    between the two centers of mass)*/
    ctmsm msmrdMarkovModel::getMSM(std::string key) {
        auto search = ctmsmDictionary.find(key);
        if (search != ctmsmDictionary.end()) {
            return search->second;
        } else {
            throw std::range_error("Key does not exist in dictionary");
        }
    }

}
