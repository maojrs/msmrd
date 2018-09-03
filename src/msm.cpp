//
// Created by maojrs on 7/4/18.
//
#include <iostream>
#include <cstdlib>
#include <ctime>
#include "particle.hpp"
#include "msm.hpp"

/**
 *  Implementations for abstract msm class inherited by all child classes
 *  @msmid Id of the Markov model to match with particle types in simulations
 *  @param tempmatrix reference to the transition rate/transition probability matrix
 *  @param lagtime lagtime used by Markov model, in the case of MSMs it will be a fixed
 *  value, for CTMSM it will change in each iteration.
 *  @param seed variable for random number generation (Note seed = -1 corresponds to random device)
 */
msmbase::msmbase(int msmid,  std::vector<std::vector<double>> &tempmatrix, double lagtime, long seed)
        : msmid(msmid), lagtime(lagtime), seed(seed){

    // Resize vectors by input matrix size and set seed of random number generator
    nstates = static_cast<int>(tempmatrix.size());
    Dlist.resize(nstates);
    Drotlist.resize(nstates);
    tmatrix.resize(nstates);
    randg.setSeed(seed);

    // Verify input 2D vector is a square matrix and fill tmatrix
    for (const auto& row : tempmatrix) {
        if (tempmatrix.size() != row.size()) {
            throw std::range_error("MSM matrix must be a square matrix");
        }
    }
    for (int i=0; i<nstates; i++) {
        tmatrix[i].resize(nstates);
        std::copy_n(tempmatrix[i].begin(), nstates, tmatrix[i].begin());
    }
};


/**
 * Implementation of discrete-time msm (msm) class, see msmbase parent class for parameter description.
 */
msm::msm(int msmid,  std::vector<std::vector<double>> &tempmatrix, double lagtime, long seed)
        : msmbase(msmid,  tempmatrix, lagtime, seed) {
    // Verify MSM transition matrix rows sum to 1 and components between 0 and 1
    double rowsum;
    for (const auto& row : tempmatrix) {
        rowsum = 0;
        for (auto& n : row){ rowsum += n; }
        if (std::abs(rowsum - 1) > tolerance) {
            throw std::range_error("Discrete-time MSM transition matrix rows should sum to 1");
        }
        for(const auto& element : row ) {
            if (element < 0 || element >1) {
                throw std::range_error("Elements of transition matrix must be probabilities (between 0 and 1)");
            }
        }
    }
};

// MISSING IMPLEMENTATION
void msm::propagate(particleMS &part, int ksteps) {

};


/**
 * Implementation of continuous-time msm (ctmsm) class, see msmbase parent class for parameter description.
 */
ctmsm::ctmsm(int msmid,  std::vector<std::vector<double>> &tempmatrix, long seed)
        : msmbase(msmid,  tempmatrix, 0.0, seed) {
    // Verify CTMSM transition rate matrix rows sum to 0
    double long rowsum;
    for (const auto& row : tempmatrix) {
        rowsum = 0;
        for (auto& n : row){ rowsum += n; }
        if (std::abs(rowsum) > tolerance) {
            throw std::range_error("Continuous-time MSM transition rate matrix rows should sum to 0");
        }
    }
    calculateParameters();
    lagtime = 0;
};

// Calculates parameters often used by ctmsm::propagate from transition matrix
void ctmsm::calculateParameters() {
    std::vector<double> ratevector(nstates-1);
    int index;
    lambda0.resize(nstates);
    ratescumsum.resize(nstates);
    // Calculates lambda0 (sum of outgoing rates), for each state
    for (int row = 0; row < nstates; row++) {
        index = 0;
        for (int col = 0; col < nstates; col++) {
            if (col != row) {
                ratevector[index] = tmatrix[row][col];
                index++;
            }
        }
        lambda0[row] = 0;
        for (auto& n : ratevector){ lambda0[row] += n; }
        // Calculates ratescumsum (rate cumulative sum for current row)
        ratescumsum[row].resize(nstates-1);
        std::copy_n(ratevector.begin(), nstates, ratescumsum[row].begin());
        for (int col = 1; col < nstates-1; col++) {
            ratescumsum[row][col] += ratescumsum[row][col-1];
        }
    }
};

/* Propagates CTMSM using the Gillespie algorithm without updating the state in the particle, useful for
 * integration with diffusion and rotation integrator */
void ctmsm::propagateNoUpdate(particleMS &part,int ksteps) {
    double lagt = 0;
    double r1, r2;
    int state = 0;
    int currentState = 1*part.state;
    for (int m = 0; m < ksteps; m++) {
        // Begins Gillespie algorithm, calculates which transition and when will it occur.
        r1 = randg.uniformRange(0,1);
        r2 = randg.uniformRange(0,1);
        lagt += std::log(1.0 / r1) / lambda0[currentState];
        for (int col = 0; col < nstates; col++) {
            if (r2 * lambda0[currentState] <= ratescumsum[currentState][col]){
                if (col < currentState) {
                    state = col;
                } else {
                    state = col + 1;
                }
                break;
            }
        }
        part.setNextState(state);
        currentState = 1*state;
        part.setLagtime(lagt);
        lagtime = lagt;
    }
};

// Propagates CTMSM using the Gillespie algorithm, updates state immediately.
void ctmsm::propagate(particleMS &part,int ksteps) {
    double lagt = 0;
    double r1, r2;
    int state = 0;
    for (int m = 0; m < ksteps; m++) {
        // Begins Gillespie algorithm, calculates which transition and when will it occur.
        r1 = randg.uniformRange(0,1);
        r2 = randg.uniformRange(0,1);
        lagt += std::log(1.0 / r1) / lambda0[part.state];
        for (int col = 0; col < nstates; col++) {
            if (r2 * lambda0[part.state] <= ratescumsum[part.state][col]){
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
};

