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
     * Trajectory classes for patchy dimer trajectories. Patchy dimer here refers to two
     * patchy particles with two patches.
     */

    /* Patchy dimer trajectory will have one unbound state (0) and eight bound states. The eight
     * bound states arise from from the four possible ways to bind through the patches combined
     * with two stable relative orientations. Used as a first toy example model to implement MSM/RD.
     * Use with patchyParticleAngular potential with numAngularMinima = 2. */
    class patchyDimer : public discreteTrajectory<8> {
    public:

        patchyDimer(unsigned long Nparticles, int bufferSize);

        patchyDimer(unsigned long Nparticles, int bufferSize, double rLowerBound, double rUpperBound);

        void setBoundStates() override;

    };


    /* Patchy dimer trajectory will have one unbound state (0) and four bound states. Same as patchyDimer
     * but with only one stable relative orientation between the particles. Used to construct
     * the base MSM model for the pentamer formation. Use with patchyParticleAngular potential
     * with numAngularMinima = 1. */
    class patchyDimer2 : public discreteTrajectory<4> {
    public:

        patchyDimer2(unsigned long Nparticles, int bufferSize);

        patchyDimer2(unsigned long Nparticles, int bufferSize, double rLowerBound, double rUpperBound);

        void setBoundStates() override;

    };


}
