//
// Created by dibakma on 27.06.18.
//
#include <array>
#include <utility>
#include "integrators/integrator.hpp"
#include "particle.hpp"

/**
 * Implementation of integrator abstract class inherited by all child classes
 */
integrator::integrator(double dt, long seed, bool rotation)
        : dt(dt), seed(seed), rotation(rotation) {
    nullPotential nullpot;
    extPotential = &nullpot; // set potential to zero as default
    randg.setSeed(seed);
    clock = 0;
};

// Function to set custom potential function
void integrator::setExternalPotential(externalPotential *pot) {
    extPotential = pot;
}

 // Integrate list of particles instead of single one (need to override in case of interacting or MS particles)
void integrator::integrateList(std::vector<particle> &parts) {
    for (int i=0; i<parts.size(); i++) {
        integrateOne(parts[i]);
    }
    clock += dt;
}



