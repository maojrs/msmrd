//
// Created by maojrs on 2/6/19.
//

#include "markovModels/continuousTimeMarkovModel.hpp"


namespace msmrd {
    /**
     * Implementation of continuous-time msm (ctmsm) class, see msmbase parent class for parameter description.
     */
    continuousTimeMarkovStateModel::continuousTimeMarkovStateModel(int msmid,
                                                                   std::vector<std::vector<double>> &tempmatrix,
                                                                   long seed)
            : markovModel(msmid, tempmatrix, 0.0, seed) {
        // Verify CTMSM transition rate matrix rows sum to 0
        double long rowsum;
        for (const auto &row : tempmatrix) {
            rowsum = 0;
            for (auto &n : row) { rowsum += n; }
            if (std::abs(rowsum) > tolerance) {
                throw std::range_error("Continuous-time MSM transition rate matrix rows should sum to 0");
            }
        }
        calculateParameters();
        lagtime = 0;
    };

    // Calculates parameters often used by ctmsm::propagate from transition matrix
    void continuousTimeMarkovStateModel::calculateParameters() {
        std::vector<double> ratevector(nstates - 1);
        int index;
        lambda0.resize(nstates);
        ratescumsum.resize(nstates);
        // Calculates lambda0 (sum of outgoing rates), for each state (sum of row excluding diagonal negative value)
        for (int row = 0; row < nstates; row++) {
            index = 0;
            for (int col = 0; col < nstates; col++) {
                if (col != row) {
                    ratevector[index] = tmatrix[row][col];
                    index++;
                }
            }
            lambda0[row] = 0;
            for (auto &n : ratevector) { lambda0[row] += n; }
            // Calculates ratescumsum (rate cumulative sum for current row)
            ratescumsum[row].resize(nstates - 1);
            std::copy_n(ratevector.begin(), nstates, ratescumsum[row].begin());
            for (int col = 1; col < nstates - 1; col++) {
                ratescumsum[row][col] += ratescumsum[row][col - 1];
            }
        }
    };

    /* Propagates CTMSM using the Gillespie algorithm without updating the state in the particle, useful for
     * integration with diffusion and rotation integrator */
    void continuousTimeMarkovStateModel::propagateNoUpdate(particleMS &part, int ksteps) {
        double lagt = 0;
        double r1, r2;
        int state = 0;
        int currentState = 1 * part.state;
        for (int m = 0; m < ksteps; m++) {
            // Begins Gillespie algorithm, calculates which transition and when will it occur.
            r1 = randg.uniformRange(0, 1);
            r2 = randg.uniformRange(0, 1);
            lagt += std::log(1.0 / r1) / lambda0[currentState];
            for (int col = 0; col < nstates; col++) {
                if (r2 * lambda0[currentState] <= ratescumsum[currentState][col]) {
                    if (col < currentState) {
                        state = col;
                    } else {
                        state = col + 1;
                    }
                    break;
                }
            }
            part.setNextState(state);
            currentState = 1 * state;
            part.setLagtime(lagt);
            lagtime = lagt;
        }
    };

    // Propagates CTMSM using the Gillespie algorithm, updates state immediately.
    void continuousTimeMarkovStateModel::propagate(particleMS &part, int ksteps) {
        double lagt = 0;
        double r1, r2;
        int state = 0;
        for (int m = 0; m < ksteps; m++) {
            // Begins Gillespie algorithm, calculates which transition and when will it occur.
            r1 = randg.uniformRange(0, 1);
            r2 = randg.uniformRange(0, 1);
            lagt += std::log(1.0 / r1) / lambda0[part.state];
            for (int col = 0; col < nstates; col++) {
                if (r2 * lambda0[part.state] <= ratescumsum[part.state][col]) {
                    if (col < part.state) {
                        state = col;
                    } else {
                        state = col + 1;
                    }
                    break;
                }
            }
            part.setState(state);
            part.setNextState(state);
            part.setLagtime(lagt);
            lagtime = lagt;
        }
    }
    
}

