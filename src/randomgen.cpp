//
// Created by maojrs on 7/25/18.
//
#include "randomgen.hpp"
#include "vec3.hpp"

/**
 * Functions of randomgen class
 */

// Sets seed for random generator mt19937_64
void randomgen::setSeed(long newseed){
    if (newseed == -1){
        mt_rand = std::mt19937_64(rd()); //random device seed
    } else {
        mt_rand = std::mt19937_64(newseed); //fixed seed
    }
};

// Returns random number between rmin and rmax sampled uniformly
double randomgen::uniformRange(double rmin, double rmax){
    std::uniform_real_distribution<double> uniform(rmin,rmax);
    return uniform(mt_rand);
};

// Returns random number sampled from normal distribution with given mean and standard deviation
double randomgen::normal(double mean, double stddev){
    std::normal_distribution<double> normaldist (mean,stddev);
    return normaldist(mt_rand);
};

// Returns random 3D vector sampled from 3D normal distribution with given mean and standard deviation
vec3<double> randomgen::normal3D(double mean, double stddev) {
    vec3<double> randvec;
    std::normal_distribution<double> normaldist (mean,stddev);
    randvec[0] = normaldist(mt_rand);
    randvec[1] = normaldist(mt_rand);
    randvec[2] = normaldist(mt_rand);
    return randvec;
};
