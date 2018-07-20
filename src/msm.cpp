//
// Created by maojrs on 7/4/18.
//
#include <functional>
#include "particle.hpp"
#include "msm.hpp"
#include <iostream>


#include <cstdlib>
#include <ctime>

/**
 * Implementation of functions for discrete-time msm (dtmsm) and
 * continuous-time msm (ctmsm) classes.
 */

int msm::propagate(particle part, int ksteps) {

};

// Calculate once the parameters often used by ctmsm::propagate from transition matrix
void ctmsm::calculateParameters() {
    std::vector<double> ratevector(nstates-1);
    int index;
    _lambda0.resize(nstates);
    _ratescumsum.resize(nstates);
    // Calculates lambda0, sum of outgoing rates, for each state
    for (int row = 0; row < nstates; row++) {
        index = 0;
        for (int col = 0; col < nstates; col++) {
            if (col != row) {
                ratevector[index] = tmatrix[row][col];
                index++;
            }
        }
        _lambda0[row] = std::accumulate(ratevector.rbegin(), ratevector.rend(), 0); //sum all elements of vector
        // Calculates rates cumulative sum for current row
        _ratescumsum[row].resize(nstates-1);
        std::copy_n(ratevector.begin(), nstates, _ratescumsum[row].begin());
        for (int col = 1; col < nstates-1; col++) {
            _ratescumsum[row][col] += _ratescumsum[row][col-1];
        }
    }
    _paramsCalculated = true;
    };

// Propagate using the Gillespie algorithm
int ctmsm::propagate(particle part,int ksteps) {
    lagtime = 0;
    double r1, r2;
    int state = 0;
    if (!_paramsCalculated) { calculateParameters(); } // Calculate parameters if not done before
    for (int m = 0; m < ksteps; m++) {
        // Begins Gillespie algorithm, calculates which transition and when will it occur.
        r1 = 1.0*std::rand()/RAND_MAX; //mt_rand();
        r2 = 1.0*std::rand()/RAND_MAX; //mt_rand();
        lagtime += std::log(1.0 / r1) / _lambda0[part.state];
        for (int col = 0; col < nstates; col++) {
            if (r2 * _lambda0[part.state] <= _ratescumsum[part.state][col]){
                if (col < part.state) {
                    state = col;
                } else {
                    state = col + 1;
                }
                break;
            }
        }
    }
    return state;
};

//template <typename T>
//ExampleThree<T>::ExampleThree(T value) : _value(value) {};
