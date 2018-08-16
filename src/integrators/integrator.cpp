//
// Created by dibakma on 27.06.18.
//
#include <array>
#include <utility>
#include "integrators/integrator.hpp"
#include "particle.hpp"
#include "msm.hpp"
#include "potentials.hpp"

/**
 * Implementations for abstract class inherited by all child classes
 */
 
 // Integrate list of particles instead of single one (need to be overriden for interacting or MS particles)
void integrator::integrateList(std::vector<particle> &parts) {
    for (int i=0; i<parts.size(); i++) {
        integrateOne(parts[i]);
    }
    clock += dt;
};



