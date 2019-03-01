//
// Created by maojrs on 1/22/19.
//

#pragma once
#include <memory>
#include "trajectory.hpp"
#include "spherePartition.hpp"
#include "tools.hpp"

namespace msmrd {
    /**
     * Trajectory class for patchy  dimer trajectories. Patchy dimer here refers to two
     * patchy particles with two patches. The main difference with the general class
     * trajectoryPositionOrientation is that this class can sample the discrete trajectory that is
     * specific for each example/application. In short words,  it chooses how to discretize the full
     * trajectory into a discretized trajectory to be analyzed and extracted into
     * a Markov state model.
     *
     * Patchy dimer will have one unbound state (0), two bound states (1,2) and
     * angularStates = numSections(numSections+1)/2 angular discrete states (3-3+angularStates).
     */
    class patchyDimer : public trajectoryPositionOrientation {
    private:

        std::vector<quaternion<double>> rotMetastableStates;
        std::vector<quaternion<double>> symmetryQuaternions{{1,0,0,0}};
        std::unique_ptr<spherePartition> spherePart;
        int angularStates;
        double tolerance = 0.25;
        /*
         * @param rotMetastableStates[X] correspond to the list of relative rotations that correspond to a
         * metastable state/region X. Each rotation is represented by a quaternion.
         * @symmetryQuaternions list of rotation that represents the structural symmetries of the patchy particle if
         * needed to define the states.
         * @param spherePart pointer to equal area spherical partition class with relevant functions
         * @param number of angularStates for discretization
         * @ param tolerance is the maximum distance away from metastable state to still be considered metastable.
         */
    public:


        patchyDimer(unsigned long Nparticles, int bufferSize, int numPartitions);

        void sampleDiscreteTrajectory(double time, std::vector<particle> &particleList) override;

        void setMetastableRegions();

        int getBoundState(quaternion<double> q1, quaternion<double> q2);

    };

}
