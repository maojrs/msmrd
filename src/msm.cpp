//
// Created by maojrs on 7/4/18.
//
#include <iostream>
#include <cstdlib>
#include <ctime>
#include "particle.hpp"
#include "msm.hpp"

/**
 * Implementation of functions for discrete-time msm (msm) and
 * continuous-time msm (ctmsm) classes.
 */

//Discrete-time msm (msm)
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

void msm::propagate(particle &part, int ksteps) {

};


//Continuous-time msm (ctmsm)
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

// Propagates CTMSM using the Gillespie algorithm
void ctmsm::propagate(particle &part,int ksteps) {
    lagtime = 0;
    double r1, r2;
    int state = 0;
    for (int m = 0; m < ksteps; m++) {
        // Begins Gillespie algorithm, calculates which transition and when will it occur.
        r1 = randg.uniformRange(0,1);
        r2 = randg.uniformRange(0,1);
        lagtime += std::log(1.0 / r1) / lambda0[part.state];
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
    }
};

