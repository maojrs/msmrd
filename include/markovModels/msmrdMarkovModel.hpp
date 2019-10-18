//
// Created by maojrs on 8/19/19.
//

#pragma once
#include "markovModels/discreteTimeMarkovModel.hpp"
#include "randomgen.hpp"

using msm = msmrd::discreteTimeMarkovStateModel;

namespace msmrd {

    class msmrdMarkovModel : public msm {
    protected:
        int numBoundStates;
        unsigned int maxNumberBoundStates;
    public:
        std::vector<int> activeSet;
        /**
         * Discrete time Markov model to be used by MSM/RD algorithm. The MSM only keeps track of the active sets,
         * so there is a difference between the indexing of the MSM and the actual state numbering.
         * @param numBoundStates total number of bound states in the two particle system
         * @param maxNumberBoundStates maximum number of possible bound states. It needs to be specified for indexing
         * purposes. After this index number, the states correspond to transition states.
         * @param activeSet corresponds to the indexes of the active states from the original discrete trajectory.
         * The transition matrix will only be as large as the number of active sets/states, so the indexing of the
         * transition matrix (tmarix) and the indexing of the discrete trajectories doesn't match. The ith component
         * of this vector yields its corresponding index in the original discrete trajectory. This is usually given
         * by pyemma when calculating an MSM.
         * @param tmatrix transition probability matrix of the MSM. Note the rank is larger than numBoundStates
         * since it also includes all the transition states
         * @param lagtime the lagtime of the discrete time MSM. This should be given in units corresponding to
         * the ones used in the simulation
         */

        msmrdMarkovModel(int numBoundStates, int maxNumberBoundStates,
                                 std::vector <std::vector<double>> tmatrix, std::vector<int> activeSet,
                                 double lagtime, long seed);

        std::tuple<double, int> calculateTransition(int initialState);

        int getMSMRDindex(int MSMindex);

        int getMSMindex(int MSMRDindex);

        void setMaxNumberBoundStates(unsigned int newIndex) {
            maxNumberBoundStates = newIndex;
        }

        unsigned int getMaxNumberBoundStates() {
            return maxNumberBoundStates;
        }

        unsigned int getNumberBoundStates() {
            return numBoundStates;
        }

        void setDbound(std::vector<double> &D, std::vector<double> &Drot);

    };

}
