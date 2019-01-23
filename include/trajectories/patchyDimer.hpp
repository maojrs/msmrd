//
// Created by maojrs on 1/22/19.
//

#pragma once
#include "trajectory.hpp"

namespace msmrd {
    /**
     * Trajectory class for patchy  dimer trajectories. Patchy dimer here refers to two
     * patchy particles with two patches. The main difference with the general class
     * trajectoryPositionOrientation is the sampling of the discrete trajectory that is
     * specific for each example/application.
     */
    class patchyDimer : public trajectoryPositionOrientation {
    public:

        // Inherit constructor from parent class
        using trajectoryPositionOrientation::trajectoryPositionOrientation;

        void sampleDiscreteTrajectory(double time, std::vector<particle> &particleList);
    };

}
