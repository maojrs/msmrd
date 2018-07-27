//
// Created by maojrs on 7/25/18.
//
#pragma once
#include <random>

#ifndef MSMRD2_RANDOMGEN_HPP
#define MSMRD2_RANDOMGEN_HPP
#endif //MSMRD2_RANDOMGEN_HPP

/**
 * Class to generate good random numbers (mt_19937_64 generator)
 * and sample from common distirbutions
 */

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
