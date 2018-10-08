//
// Created by dibakma on 22.06.18.
//
#pragma once

#include <array>
#include "particle.hpp"
#include "integrators/integrator.hpp"
#include "trajectories/trajectory.hpp"

namespace msmrd {
    class simulation {
    public:
        integrator& integ;
        std::vector<particle> &particleList;
        simulation(integrator& integ, std::vector<particle> &particleList) : integ(integ), particleList(particleList) {};
        void run(const int steps, trajectory& traj, int stride);
    };

}