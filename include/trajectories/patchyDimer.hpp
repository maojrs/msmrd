//
// Created by maojrs on 1/22/19.
//

#pragma once
#include "trajectory.hpp"
#include "tools.hpp"

namespace msmrd {
    /**
     * Trajectory class for patchy  dimer trajectories. Patchy dimer here refers to two
     * patchy particles with two patches. The main difference with the general class
     * trajectoryPositionOrientation is the sampling of the discrete trajectory that is
     * specific for each example/application.
     */
    class patchyDimer : public trajectoryPositionOrientation {
    private:
        /*
         * @param rotMetastableStates[X] correspond to the list of equivalent rotations that yield a certain
         * metastable state/region X. Each rotation is represented by a quaternion
         */
        std::vector<std::vector<quaternion<double>>> rotMetastableStates;

    public:

        // Inherit constructor from parent class
        using trajectoryPositionOrientation::trajectoryPositionOrientation;

        void sampleDiscreteTrajectory(double time, std::vector<particle> &particleList);

        void setMetastableRegions();
    };

}
