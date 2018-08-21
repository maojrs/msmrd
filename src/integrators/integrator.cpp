//
// Created by dibakma on 27.06.18.
//
#include <array>
#include <utility>
#include "integrators/integrator.hpp"
#include "particle.hpp"

/**
 * Implementations for abstract class inherited by all child classes
 */

// Constructor that stes seed and potential to null
integrator::integrator(double dt, long seed, bool rotation)
        : dt(dt), seed(seed), rotation(rotation) {
    nullPotential nullpot;
    extPotential = &nullpot;
    randg.setSeed(seed);
    clock = 0;
};

// fucntion to set potential
void integrator::setExternalPotential(externalPotential *pot) {
    extPotential = pot;
}

 // Integrate list of particles instead of single one (need to be overriden for interacting or MS particles)
void integrator::integrateList(std::vector<particle> &parts) {
    for (int i=0; i<parts.size(); i++) {
        integrateOne(parts[i]);
    }
    clock += dt;
}



