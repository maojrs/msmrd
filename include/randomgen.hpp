//
// Created by maojrs on 7/25/18.
//
#pragma once
#include <random>

#ifndef MSMRD2_RANDOMGEN_HPP
#define MSMRD2_RANDOMGEN_HPP
#endif //MSMRD2_RANDOMGEN_HPP

class randomgen {
private:
    long seed;
    std::random_device rd;
    std::mt19937_64 mt_rand;  // random number generator
public:
    randomgen() = default;

    void setSeed(long newseed);
    double uniformRange(double rmin, double rmax);
    double normal(double mean, double stddev);
};
