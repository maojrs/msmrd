//
// Created by maojrs on 7/4/18.
//
#include <iostream>
#include <cstdlib>
#include <ctime>
#include "particle.hpp"
#include "markovModels/markovModel.hpp"


namespace msmrd {
    /**
     *  Implementations for abstract msm class inherited by all child classes
     *  @msmid Id of the Markov model to match with particle types in simulations
     *  @param tempmatrix reference to the transition rate/transition probability matrix
     *  @param lagtime lagtime used by Markov model, in the case of MSMs it will be a fixed
     *  value, for CTMSM it will change in each iteration.
     *  @param seed variable for random number generation (Note seed = -1 corresponds to random device)
     */
    markovModel::markovModel(int msmid, std::vector<std::vector<double>> &tempmatrix, double lagtime, long seed)
            : msmid(msmid), lagtime(lagtime), seed(seed) {

        // Resize vectors by input matrix size and set seed of random number generator
        nstates = static_cast<int>(tempmatrix.size());
        Dlist.resize(nstates);
        Drotlist.resize(nstates);
        tmatrix.resize(nstates);
        randg.setSeed(seed);

        // Verify input 2D vector is a square matrix and fill tmatrix
        for (const auto &row : tempmatrix) {
            if (tempmatrix.size() != row.size()) {
                throw std::range_error("MSM matrix must be a square matrix");
            }
        }
        for (int i = 0; i < nstates; i++) {
            tmatrix[i].resize(nstates);
            std::copy_n(tempmatrix[i].begin(), nstates, tmatrix[i].begin());
        }
    };

}
