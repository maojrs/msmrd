//
// Created by maojrs on 2/6/19.
//

#pragma once
#include <map>
#include "markovModels/continuousTimeMarkovModel.hpp"

namespace msmrd {
    using ctmsm = msmrd::continuousTimeMarkovStateModel;
    /*
     * Base Markov state model class used for MSM/RD algorithm. Although, it is based on the continuous-time
     * Markov state model class, its functionality and construction differs enough such that it is not convenient to
     * inherit from the main class. Actually the Markov model is completely constructed with the rate dictionary, and
     * it defines a large list/dictionary of ctmsms, one for each discrete region in order to emulate the spatially
     * dependent rates.
     */
    class msmrdMarkovModel {
    protected:
        long seed;
        std::map<std::string, float> rateDictionary;
        std::map<std::string, ctmsm> ctmsmDictionary = {};
    public:
        unsigned int nstates;
        unsigned int numDiscreteOrientations;
        /**
        * @param seed variable for random number generation (Note values of seed <= -1 correspond to random device)
        * @param rateDictionary dictionary relating transitions to its corresponding rates. The keys
        * are strings of the form "state1->state2". The bound states have "b" before the number, and the transition
        * states are denoted by integers.
        * @param ctmsmDictionary dicitionary mapping keys of relevant transitions to the corresponding
        * continuous time Markov model.
        * @param nstates number of states in the msm (including unbound state)
        * @param numDiscreteOrientations number of Galerkin discrete orientations of one particle.
        */

        msmrdMarkovModel(unsigned int nstates, unsigned int numDiscreteOrientations,
                         long seed, std::map<std::string, float> &rateDictionary);

        void generateMarkovModels();

        float getRate(std::string key);

        ctmsm getMSM(std::string key);

    };

}


