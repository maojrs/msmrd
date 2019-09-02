//
// Created by maojrs on 1/22/19.
//

#pragma once
#include <memory>
#include "trajectories/discrete/discreteTrajectory.hpp"
#include "discretizations/positionOrientationPartition.hpp"
#include "tools.hpp"
#include "H5Cpp.h"

namespace msmrd {
    /**
     * Trajectory class for patchy dimer trajectories. Patchy dimer here refers to two
     * patchy particles with two patches. The main difference with the general class
     * trajectoryPositionOrientation is that this class can sample the discrete trajectory that is
     * specific for the patchy dimer application. In short words,  it chooses how to discretize the full
     * trajectory of two particle into a discretized trajectory to be analyzed and extracted into
     * a Markov state model. It also implements functionality to discretize trajectories directly loaded
     * from a python array.
     *
     * This Patchy dimer trajectory will have one unbound state (0), eight bound states.
     */
    class patchyDimer : public discreteTrajectory<8> {
    public:

        patchyDimer(unsigned long Nparticles, int bufferSize);

        patchyDimer(unsigned long Nparticles, int bufferSize, double rLowerBound, double rUpperBound);

        void setBoundStates() override;


    };

}
