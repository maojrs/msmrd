//
// Created by maojrs on 2/6/19.
//

#include "markovModels/continuousTimeMarkovModel.hpp"


namespace msmrd {
    /**
     * Implementation of continuous-time msm (ctmsm) class, see msmbase parent class for parameter description.
     */
    continuousTimeMarkovStateModel::continuousTimeMarkovStateModel(int msmid, long seed)
        : markovModel(msmid, 0.0, seed) {
        calculateParameters();
        lagtime = 0;
    };

    continuousTimeMarkovStateModel::continuousTimeMarkovStateModel(int msmid,
                                                                   std::vector<std::vector<double>> tmatrix,
                                                                   long seed)
            : markovModel(msmid, tmatrix, 0.0, seed) {
        // Verify CTMSM transition rate matrix rows sum to 0
        double long rowsum;
        for (const auto &row : tmatrix) {
            rowsum = 0;
            for (auto &n : row) { rowsum += n; }
                if (std::abs(rowsum) > tolerance) {
                    throw std::invalid_argument("Continuous-time MSM transition rate matrix rows should sum to 0");
                }
        }
        calculateParameters();
        lagtime = 0;
    };

    // Calculates parameters often used by ctmsm::propagate from transition matrix
    void continuousTimeMarkovStateModel::calculateParameters() {
        lambda0.resize(nstates);
        ratescumsum.resize(nstates, std::vector<double>(nstates - 1));
        if (nstates == 1) {
            lambda0[0] = 0;
            ratescumsum[0] = {};
        }
        else { //nstates > 1
            std::vector<double> ratevector;
            ratevector.resize(nstates - 1);
            int index;
            //Calculates lambda0 (sum of outgoing rates), for each state (sum of row excluding diagonal negative value)
            for (int row = 0; row < nstates; row++) {
                index = 0;
                // ratevector will be the current row vector without the diagonal component
                for (int col = 0; col < nstates; col++) {
                    if (col != row) {
                        ratevector[index] = tmatrix[row][col];
                        index++;
                    }
                }
                // lambda 0  is the sum of all components of ratevector
                lambda0[row] = std::accumulate(ratevector.begin(), ratevector.end(),0);
                // Calculates ratescumsum (rate cumulative sum for current row)
                std::copy_n(ratevector.begin(), nstates - 1, ratescumsum[row].begin());
                for (int col = 1; col < nstates - 1; col++) {
                    ratescumsum[row][col] += ratescumsum[row][col - 1];
                }
            }
        }
    };

    /* Propagates CTMSM using the Gillespie algorithm without updating the state in the particle, useful for
     * integration with diffusion and rotation integrator */
    void continuousTimeMarkovStateModel::propagateNoUpdate(particle &part, int ksteps) {
        if (nstates > 1) {
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
        }
        else {
            throw std::invalid_argument("Nothing to propagate since the rate matrix is a scalar (0)");
        }
    }

    // Propagates CTMSM in particles using the Gillespie algorithm, updates state  of particles immediately.
    void continuousTimeMarkovStateModel::propagate(particle &part, int ksteps) {
        if (nstates > 1) {
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
        else {
            throw std::invalid_argument("Nothing to propagate since the rate matrix is a scalar (0)");
        }
    }

    /* Propagates CTMSM given an initial state for ksteps timesteps. Returns tuple with total lagtime
     * and final state. Unlike propagate fucntion, this function doesn't involve any particles.*/
    std::tuple<double, int> continuousTimeMarkovStateModel::propagateMSM(int initialState, int ksteps) {
        double lagt = 0;
        double r1, r2;
        int state = initialState;
        int newState;
        if (nstates > 1) {
            for (int m = 0; m < ksteps; m++) {
                // Begins Gillespie algorithm, calculates which transition and when will it occur.
                r1 = randg.uniformRange(0, 1);
                r2 = randg.uniformRange(0, 1);
                lagt += std::log(1.0 / r1) / lambda0[state];
                for (int col = 0; col < nstates; col++) {
                    if (r2 * lambda0[state] <= ratescumsum[state][col]) {
                        if (col < state) {
                            newState = col;
                        } else {
                            newState = col + 1;
                        }
                        break;
                    }
                }
                state = newState;
            }
            return std::make_tuple(lagt, state);
        } else {
            throw std::invalid_argument("Nothing to propagate since the rate matrix is a scalar (0)");
        }
    }


}

