//
// Created by maojrs on 7/4/18.
//

#pragma once
#include <random>
#include <array>
#include <algorithm>

/**
 * Abstract base class for Markov state models
 */
template<unsigned int N>
class msmbase {
public:
    int msmid;
    int nstates = N;
    double lagtime;
    std::array<std::array<double, N>, N> tmatrix;
    //std::array<double, N * N> tmatrix;

    // Main constructor
    msmbase(int msmid,  std::array<std::array<double, N>, N> tmatrix, double lagtime)
            : msmid(msmid), tmatrix(tmatrix), lagtime(lagtime) {};

    // Constructor capable of receiving python numpy arrays as input
    msmbase(int msmid, std::vector<std::vector<double>> &tempmatrix, double lagtime)
            : msmid(msmid), lagtime(lagtime) {
        for (int i=0; i<N; i++) {
            std::copy_n(tempmatrix[i].begin(), N, tmatrix[i].begin());
        }
    };

    // Main functions definitions (=0 for abtract class)
    virtual void propagate() = 0;


    /** Get properties functions for pybinding **/
    int getID() { return  msmid; }
    int getNstates() { return  nstates; }
    double getLagtime() {return lagtime; }
    std::array<std::array<double, N>, N> getTmatrix() { return  tmatrix; }
    //std::array<double, N * N> getTmatrix() { return  tmatrix; }
};


/**
 * Child classes of msmbase, discrete-time (msm) and continuous-time (ctmsm)
 */

template<unsigned int N>
class msm: public msmbase<N> {
public:
    using msmbase<N>::msmbase; // constructors inheritance
    virtual void propagate() override;

};

template<unsigned int N>
class ctmsm: public msmbase<N> {
public:
    using msmbase<N>::msmbase; // constructors inheritance
};






