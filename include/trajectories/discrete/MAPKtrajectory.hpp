//
// Created by maojrs on 5/3/21.
//

#pragma once
#include <memory>
#include "trajectories/discrete/discreteTrajectory.hpp"
#include "discretizations/positionOrientationPartition.hpp"
#include "tools.hpp"

namespace msmrd {
/**
 * Trajectory class for patchy protein trajectory for the MAPK model (specific application). This class
 * samples the discrete trajectory. We just need to specify how to discretize the full trajectory by
 * setting all the possible bound states between all the particles.
 *
 * Note discretization used here should match the discretization of MSM/RD for consistent results.
 */

    /* Patchy protein MAPK trajectory implementation with one unbound state (0)
     * and 4 bound states (1,2,3,4). The bound states 1 and 2 correspond to the kinase binding to sites
     * one and two, while bound states 3 and 4 correspond to the phosphatase binding to the same
     * two binding sites. */
    class MAPKtrajectory : public discreteTrajectory<4> {
    public:

        MAPKtrajectory(unsigned long Nparticles, int bufferSize);

        MAPKtrajectory(unsigned long Nparticles, int bufferSize, double rLowerBound, double rUpperBound);

        void setBoundStates();

        int sampleDiscreteState(particle part1, particle part2) override;

    };


}