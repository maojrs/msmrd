//
// Created by maojrs on 7/4/18.
//

#include "msm.hpp"

template<unsigned int N>
msm<N>::msm(int type, int nstates, std::array<double, N * N> tmatrix, double lagtime)
        : type(type), nstates(nstates), tmatrix(tmatrix), lagtime(lagtime) {}

template<unsigned int N>
ctmsm<N>::ctmsm(int type, int nstates, std::array<double, N * N> ratematrix, double lagtime)
        : type(type), nstates(nstates), ratematrix(ratematrix), lagtime(lagtime) {};

//void msm::propagate(){
//};