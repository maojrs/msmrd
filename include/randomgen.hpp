//
// Created by maojrs on 7/25/18.
//
#pragma once
#include <random>
#include "vec3.hpp"


/**
 * Declares class to generate good random numbers (mt_19937_64 generator)
 * and sample from common distirbutions
 */
class randomgen {
private:
    long seed = -1;
    std::random_device rd;
    std::mt19937_64 mt_rand;  // random number generator
public:
    randomgen() { mt_rand = std::mt19937_64(rd()); };
    randomgen(const randomgen&) { mt_rand = std::mt19937_64(rd()); };

    void setSeed(long newseed);
    double uniformRange(double rmin, double rmax);
    double normal(double mean, double stddev);
    vec3<double> normal3D(double mean, double stddev);
    vec3<double> uniformSphere(double maxrad);
    vec3<double> uniformShell(double minrad, double maxrad);
};
