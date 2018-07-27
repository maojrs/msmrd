//
// Created by maojrs on 7/25/18.
//
#include "randomgen.hpp"

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

// Returns random number from normal distribution with given mean and standard deviation
double randomgen::normal(double mean, double stddev){
    std::normal_distribution<double> normaldist (mean,stddev);
    return normaldist(mt_rand);
};