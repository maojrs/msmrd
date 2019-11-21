//
// Created by maojrs on 8/19/19.
//

#include "markovModels/msmrdMarkovModel.hpp"

namespace msmrd {
    /**
     * Implementation of Markov model class specializes for MSM/RD.
     */
    msmrdMarkovModel::msmrdMarkovModel(int numBoundStates, int maxNumberBoundStates,
                                                       std::vector<std::vector<double>> tmatrix,
                                                       std::vector<int> activeSet, double lagtime, long seed)
            : numBoundStates(numBoundStates), maxNumberBoundStates(maxNumberBoundStates), activeSet(activeSet),
              msm(-1, tmatrix, lagtime, seed){
        randg.setSeed(seed);
        Dlist.resize(numBoundStates);
        Drotlist.resize(numBoundStates);
    };


    /* Calculates next transition (time and next state), given an initial state using the msmrdMSM.
     * Note it returns the state in the active set indexing (see activeSet class header)*/
    std::tuple<double, int> msmrdMarkovModel::calculateTransition(int initialState) {
        // Find index of initialState in the transition matrix
        int localState = getMSMindex(initialState);
        // Return empty event if initialState is not in activeSet
        if (localState == -1) {
            return std::make_tuple(std::numeric_limits<double>::infinity(), -1);
        }
        // Calculate MSM transition
        auto transition = propagateMSM(localState, 1);
        double transitionTime = std::get<0>(transition);
        int nextState = std::get<1>(transition);
        // Recover correct indexing for MSM/RD (same as discrete trajectories)
        auto nState = activeSet[nextState];
        return std::make_tuple(transitionTime, nState);
    };


    /* This function sets the diffusion coefficients for the bound states of A and B.
    * Note size of vectors must match numBoundStates. This functions NEEDS
    * to be called to fill in the diffusion coefficients. */
    void msmrdMarkovModel::setDbound(std::vector<double> &D, std::vector<double> &Drot) {
        if ( (D.size() != numBoundStates ) or (Drot.size() != numBoundStates) ) {
            throw std::invalid_argument("Vectors of diffusion coefficients must match number of bound states");
        }
        Dlist = D;
        Drotlist = Drot;
    }


    //  Get MSMRD state from index
    int msmrdMarkovModel::getActiveSetIndex(int MSMindex){
        return activeSet[MSMindex];
    };

    int msmrdMarkovModel::getMSMindex(int activeSetIndex){
        // Find index (MSMindex) of initialState (MSMRDindex) in the transition matrix
        std::vector<int>::iterator itr = std::find(activeSet.begin(), activeSet.end(), activeSetIndex);
        if (itr != activeSet.cend()) {
            // Element found, return corresponding MSM index
            return std::distance(activeSet.begin(), itr);
        }
        else {
            // Element not found, return -1
            return -1;
        }
    }


}
