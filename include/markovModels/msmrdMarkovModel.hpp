//
// Created by maojrs on 2/6/19.
//

#pragma once
#include <map>
#include "markovModels/continuousTimeMarkovModel.hpp"

namespace msmrd {
    using ctmsm = msmrd::continuousTimeMarkovStateModel;
    /*
     * Base Markov state model class used for MSM/RD algorithm. This Markov model is completely constructed from the
     * rate dictionary. By default the rate dictionary passes only the rates from transition states to bound states.
     * However, if all the rates between all the states is specified, it can further define a large list/dictionary
     * of ctmsms, one for each discrete region in order to emulate the spatially dependent rates. However, this latter
     * feature is optional for alternative formulations of the scheme. In general the rate dictionary of transitions
     * between transition and bound states will be enough.
     */
    class msmrdMarkovModel {
    protected:
        long seed;
        std::map<std::string, float> rateDictionary;
        randomgen randg;
    public:
        unsigned int numBoundStates;
        unsigned int numTransitionStates;
        //unsigned int numDiscreteOrientations;
        /**
        * @param seed variable for random number generation (Note values of seed <= -1 correspond to random device)
        * @param rateDictionary dictionary relating transitions to its corresponding rates. The keys
        * are strings of the form "state1->state2". The bound states have "b" before the number, and the transition
        * states are denoted by integers.
        * @param random number generator class based on mt19937
        * @param numBoundStates number of bound states in the msm
        * @param numTransitionStates number of transition states
        */

        msmrdMarkovModel(unsigned int numBoundStates, unsigned int numTransitionStates,
                         long seed, std::map<std::string, float> &rateDictionary);

        std::tuple<double, int> computeTransition(int transitionState);

        float getRate(std::string key);

        ctmsm getMSM(std::string key);

    };

}


