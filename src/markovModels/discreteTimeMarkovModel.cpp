//
// Created by maojrs on 2/6/19.
//

#include "markovModels/discreteTimeMarkovModel.hpp"

namespace msmrd{
    /**
    * Implementation of discrete-time msm (msm) class, see msmbase parent class for parameter description.
    */
    discreteTimeMarkovStateModel::discreteTimeMarkovStateModel(int msmid, double lagtime, long seed) :
            markovModel(msmid, 0.0, seed) {
        tmatrix = {{1.0}};
    };

    discreteTimeMarkovStateModel::discreteTimeMarkovStateModel(int msmid, std::vector<std::vector<double>> tmatrix,
                                                               double lagtime, long seed)
            : markovModel(msmid, tmatrix, lagtime, seed) {
        // Verify MSM transition matrix rows sum to 1 and components between 0 and 1
        double rowsum;
        for (const auto &row : tmatrix) {
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

    // Propagates the discrete Markov chain for ksteps
    void discreteTimeMarkovStateModel::propagate(particleMS &part, int ksteps) {
        double r1;
        double sum;
        int currentState = 1 * part.state;
        for (int m = 0; m < ksteps; m++){
            sum = 0;
            r1 = randg.uniformRange(0, 1);
            for (int i = 0; i < nstates; i++ ){
                sum += tmatrix[currentState][i];
                if (r1 <= sum) {
                    part.setState(i);
                    part.setNextState(i);
                    currentState = i;
                    break;
                }
            }
        }
    };

    /* Propagates the discrete Markov chain for ksteps, without any particles involved, returns
     * total lagtime and final state */
    std::tuple<double, int> discreteTimeMarkovStateModel::propagateMSM(int initialState, int ksteps) {
        double r1;
        double sum;
        int state = initialState;
        for (int m = 0; m < ksteps; m++){
            sum = 0;
            r1 = randg.uniformRange(0, 1);
            for (int i = 0; i < nstates; i++ ){
                sum += tmatrix[state][i];
                if (r1 <= sum) {
                    state = i;
                    break;
                }
            }
        }
        return std::make_tuple(lagtime*ksteps, state);
    };

    // Propagates the discrete Markov chain for ksteps without updating
    void discreteTimeMarkovStateModel::propagateNoUpdate(particleMS &part, int ksteps) {
        double r1;
        double sum;
        int currentState = 1 * part.state;
        for (int m = 0; m < ksteps; m++){
            sum = 0;
            r1 = randg.uniformRange(0, 1);
            for (int i = 0; i < nstates; i++ ){
                sum += tmatrix[currentState][i];
                if (sum <= r1) {
                    part.setNextState(i);
                    currentState = i;
                    break;
                }
            }
        }
    };

}