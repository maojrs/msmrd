//
// Created by dibakma on 22.06.18.
//
#pragma once

#include <array>
#include "particle.h"
#include "integrator.h"

class simulation {
public:
    std::vector<particle<double>> &particles;
    const int Nparticles;
    simulation(std::vector<particle<double>> particles): particles(particles), Nparticles(particles.size()) {};
    void run(const double timestep, const int steps);
};