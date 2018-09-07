//
// Created by dibakma on 22.06.18.
//
#pragma once

#include <array>
#include "particle.hpp"
#include "integrators/integrator.hpp"

namespace msmrd {
    class simulation {
    public:
        std::vector<particle> &particles;
        const int Nparticles;

        simulation(std::vector<particle> particles) : particles(particles), Nparticles(particles.size()) {};

        void run(const double timestep, const int steps);
    };

}