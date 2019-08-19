//
// Created by maojrs on 8/19/19.
//

#include "markovModels/msmrdMarkovModelDiscrete.hpp"

namespace msmrd {
    /**
     * Implementation of Markov model class specializes for MSM/RD.
     */
    msmrdMarkovModelDiscrete::msmrdMarkovModelDiscrete(int numBoundStates, int maxNumberBoundStates,
                                                       std::vector<std::vector<double>> tmatrix,
                                                       std::vector<int> activeSet, double lagtime, long seed)
            : numBoundStates(numBoundStates), maxNumberBoundStates(maxNumberBoundStates), activeSet(activeSet),
              msm(-1, tmatrix, lagtime, seed){
        randg.setSeed(seed);
        Dlist.resize(numBoundStates);
        Drotlist.resize(numBoundStates);
    };


    std::tuple<double, int> msmrdMarkovModelDiscrete::calculateTransition(int initialState) {
        // Find index of initialState in the transition matrix
        int localState;
        std::vector<int>::iterator itr = std::find(activeSet.begin(), activeSet.end(), initialState);
        if (itr != activeSet.cend()) {
            // Element found, local index set
            localState = std::distance(activeSet.begin(), itr);
        }
        else {
            // Element not found, so return empty transition
            return std::make_tuple(std::numeric_limits<double>::infinity(), -1);
        }
        // Calculate MSM transition
        auto transition = propagateMSM(localState, 1);
        double transitionTime = std::get<0>(transition);
        int nextState = std::get<1>(transition);
        // Recover correct indexing for MSM/RD (same as discrete trajectories)
        nextState = activeSet[nextState];
        return std::make_tuple(transitionTime, nextState);
    };


    /* This function sets the diffusion coefficients for the bound states of A and B.
    * Note size of vectors must match numBoundStates. This functions NEEDS
    * to be called to fill in the diffusion coefficients. */
    void msmrdMarkovModelDiscrete::setDbound(std::vector<double> &D, std::vector<double> &Drot) {
        if ( (D.size() != numBoundStates ) or (Drot.size() != numBoundStates) ) {
            std::__throw_range_error("Vectors of diffusion coefficients must match number of bound states");
        }
        Dlist = D;
        Drotlist = Drot;
    }


}
