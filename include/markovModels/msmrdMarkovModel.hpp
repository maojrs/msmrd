//
// Created by maojrs on 2/6/19.
//

#pragma once
#include <map>
#include "randomgen.hpp"


namespace msmrd {
    /*
     * Base Markov state model class used for MSM/RD algorithm. This Markov model is completely constructed from the
     * rate dictionary. By default the rate dictionary passes only the rates from transition states to bound states and
     * from bound states to transition states. However, if all the rates between all the states is specified, it could
     * further define a large list/dictionary of ctmsms, one for each discrete region in order to emulate the spatially
     * dependent rates. However, this latter feature is optional for alternative formulations of the scheme. In general
     * the rate dictionary of transitions between transition and bound states (and viceversa) will be enough.
     *
     * It is used to establish the coupling between unbound states and bound states and to model the transition
     * between bound states. It does NOT manage the Markov dynamics in the unbound states. That is done by the MSMs
     * in MSMlist in the msmrdIntegrator.
     */
    class msmrdMarkovStateModel {
    protected:
        long seed;
        std::map<std::string, float> rateDictionary;
        randomgen randg;
        unsigned int maxNumberBoundStates = 10;
    public:
        unsigned int numBoundStates;
        unsigned int numTransitionStates;
        std::vector<double> Dboundlist;
        std::vector<double> DboundRotlist;
        /**
        * @param seed variable for random number generation (Note values of seed <= -1 correspond to random device)
        * @param rateDictionary dictionary relating transitions to its corresponding rates. The keys
        * are strings of the form "state1->state2". The bound states have "b" before the number, and the transition
        * states are denoted by integers.
        * @param randg number generator class based on mt19937
        * @param maxNumberBoundStates maximum number of bound states supported. It is used to determine how to
        * count (index) the transition states. The state maxNumberBoundStates + 1 will correspond not to a bound state
        * but to the first transition state. This parameter has to be consistent with the one used
        * by the msmrd integrator and the with the one used to generate the state numbering in the discrete
        * trajectory (see patchyDimer trajectory for example).
        * @param numBoundStates number of bound states in the msm
        * @param numTransitionStates number of transition states
        * @param Dboundlist list of diffusion coefficients for each bound state.
        * @param DboundRotlist list of rotational diffusion coefficients for each bound state (particle type C)

        */

        msmrdMarkovStateModel(unsigned int numBoundStates, unsigned int numTransitionStates,
                         long seed, std::map<std::string, float> &rateDictionary);

        std::tuple<double, int> computeTransition2BoundState(int transitionState);

        std::tuple<double, int> computeTransition2UnboundState(int boundState);


        float getRate(std::string key);

        void setmaxNumberBoundStates(unsigned int newIndex) {
            maxNumberBoundStates = newIndex;
        }

        unsigned int getStartIndexTransitionRates() {
            return maxNumberBoundStates;
        }


        void setDbound(std::vector<double> &D, std::vector<double> &Drot);


    };

}


