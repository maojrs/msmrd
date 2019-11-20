//
// Created by maojrs on 4/8/19.
//

#pragma once
#include <memory>
#include "trajectories/discrete/discreteTrajectory.hpp"
#include "discretizations/positionOrientationPartition.hpp"
#include "tools.hpp"

namespace msmrd {
    /**
     * Trajectory classes for patchy protein trajectories. Patchy protein here refers to two
     * patchy particles each with a different number of patches. Patches can be of different types,
     * and molecules conformations can be modeled by activating/deactivating patches. These classes sample the
     * discrete trajectory, each one specific for a given patchy protein application. They just need to specify
     * how to discretize the full trajectory by setting the bound states of the two particles.
     *
     * Note sicretization used here should match the discretization of MSM/RD for consistent results.
     */

     /* Patchy protein implementation with one unbound state (0) and 6 bound states (1,2,3,4,5,6) */
    class patchyProteinTrajectory : public discreteTrajectory<6> {
    public:

        patchyProteinTrajectory(unsigned long Nparticles, int bufferSize);

        patchyProteinTrajectory(unsigned long Nparticles, int bufferSize, double rLowerBound, double rUpperBound);

        int sampleDiscreteState(particle part1, particle part2) override;

        void setBoundStates();

    };


}
