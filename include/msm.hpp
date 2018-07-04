//
// Created by maojrs on 7/4/18.
//

#pragma once
#include <random>
#include <array>

/**
 * Base classes for discrete-time and continuous-time Markov state models
 */

template<unsigned int N>
class msm {
public:
    int type;
    int nstates = N;
    double lagtime;
    std::array<double, N * N> tmatrix;

    msm(int type, int nstates, std::array<double, N * N> tmatrix, double lagtime);

    /** Get properties functions for pybinding **/
    int getType() { return  type; }
    int getNstates() { return  nstates; }
    double getLagtime() {return lagtime; }
    std::array<double, N * N> getTmatrix() { return  tmatrix; }
};

template<unsigned int N>
class ctmsm {
public:
    int type;
    int nstates = N;
    double lagtime;
    std::array<double, N * N> ratematrix;

    ctmsm(int type, int nstates, std::array<double, N * N> tmatrix, double lagtime);

    /** Get properties functions for pybinding **/
    int getType() { return  type; }
    int getNstates() { return  nstates; }
    double getLagtime() {return lagtime; }
    std::array<double, N * N> getRatematrix() { return  ratematrix; }
};

