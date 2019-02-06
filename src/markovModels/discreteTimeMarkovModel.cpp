//
// Created by maojrs on 2/6/19.
//

#include "markovModels/discreteTimeMarkovModel.hpp"

namespace msmrd{
    /**
    * Implementation of discrete-time msm (msm) class, see msmbase parent class for parameter description.
    */
    discreteTimeMarkovStateModel::discreteTimeMarkovStateModel(int msmid, std::vector<std::vector<double>> &tempmatrix, double lagtime, long seed)
            : markovModel(msmid, tempmatrix, lagtime, seed) {
        // Verify MSM transition matrix rows sum to 1 and components between 0 and 1
        double rowsum;
        for (const auto &row : tempmatrix) {
            rowsum = 0;
            for (auto &n : row) { rowsum += n; }
            if (std::abs(rowsum - 1) > tolerance) {
                throw std::range_error("Discrete-time MSM transition matrix rows should sum to 1");
            }
            for (const auto &element : row) {
                if (element < 0 || element > 1) {
                    throw std::range_error("Elements of transition matrix must be probabilities (between 0 and 1)");
                }
            }
        }
    };

    // MISSING IMPLEMENTATION
    void discreteTimeMarkovStateModel::propagate(particleMS &part, int ksteps) {

    };
    
}