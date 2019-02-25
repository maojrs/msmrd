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
     * it defines a large list/dicitionary of ctmsms, one for each discrete region in order to emulate the spatially
     * dependent rates.
     */
    class msmrdMarkovModel {
    protected:
        const long double tolerance = 1 * pow(10.0L, -10);
        long seed;
        randomgen randg;
        std::map<std::string, float> rateDictionary;
        std::map<std::string, ctmsm> ctmsmDictionary;
        //std::vector<ctmsm> ctmsmList = {};
    public:
        unsigned int nstates;
        unsigned int numDiscreteOrientations;
        double lagtime = 0;
        /**
        * @param msmid ID of the msm, corresponds to the particle type
        * @param tolerance tolerance limit for MSM integrity check
        * @param seed variable for random number generation; seed = -1 corresponds to random_device;
        * @param randg random number generator class based on mt19937
        * @param rateDictionary dictionary relating transitions to its corresponding rates. The keys
        * are strings of the form "state1->state2". The bound states hava b before the number, and the transition
        * states are denoted with duple integers.
        * @param ctmsmDictionary dicitionary mapping keys of relevant transitions to the corresponding
        * Markov model.
        * @param nstates number of states in the msm (including unbound state)
        * @param numDiscreteOrientations number of Galerkin discrete orientations of one particle.
        * @param lagtime msm lagtime (in ctmsm it is calculated after each propagation step)
        */

        msmrdMarkovModel(unsigned int nstates, unsigned int numDiscreteOrientations,
                         long seed, std::map<std::string, float> &rateDictionary);

        void generateMarkovModels();

    };

}


