//
// Created by dibakma on 27.06.18.
//

#pragma once
#include <random>
#include <array>
#include "particle.hpp"


/**
 * Base class for integrators
 */
class integrator {
protected:
    std::vector<particle>& particles;
    integrator();
    void integrate(double);
};

///**
// * Standard Brownian motion
// */
//class brownianDynamics : integrator {
//    const int nparticles;
//    std::mt19937 generator;
//    brownianDynamics(std::vector<particle>&, int);
//    void integrate(double);
//};