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
        integrator &integ;
        std::vector<particle> &particleList;
        /**
         * @param integ Integrator to be used for simulation, works for any integrator since they are all
         * childs from abstract class
         * @param particleList List of particles to be integrated. In principle it can also take lists of
         * cutom particle types as long as they are childs of the original particle class. However, watch out
         * in the bindings, since custom types might not work.
         */

        simulation(integrator &integ, std::vector<particle> &particleList);

        void run(const int steps, trajectory& traj, int stride);
    };

}