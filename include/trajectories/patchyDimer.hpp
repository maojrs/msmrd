//
// Created by maojrs on 1/22/19.
//

#pragma once
#include <memory>
#include "trajectoryPositionOrientation.hpp"
#include "discretizations/positionOrientationPartition.hpp"
#include "tools.hpp"

namespace msmrd {
    /**
     * Trajectory class for patchy  dimer trajectories. Patchy dimer here refers to two
     * patchy particles with two patches. The main difference with the general class
     * trajectoryPositionOrientation is that this class can sample the discrete trajectory that is
     * specific for the patchy dimer application. In short words,  it chooses how to discretize the full
     * trajectory of two particle into a discretized trajectory to be analyzed and extracted into
     * a Markov state model.
     *
     * Patchy dimer will have one unbound state (0), two (or more) bound states (1,2) and
     * angularStates = numSections(numSections+1)/2 angular discrete states (3-3+angularStates).
     * Note in the current indexing implementation only up to 9 bound states are supported.
     */
    class patchyDimer : public trajectoryPositionOrientation {
    private:
        std::vector<quaternion<double>> rotMetastableStates;
        std::vector<quaternion<double>> symmetryQuaternions{{1,0,0,0}};
        std::unique_ptr<spherePartition> spherePart;
        std::unique_ptr<quaternionPartition> quaternionPart;
        std::unique_ptr<positionOrientationPartition> positionOrientationPart;
        int startIndexTransitionStates = 10;
        int angularStates;
        double tolerance = 0.1;
        int prevsample = 0;
    public:
        /*
         * @param rotMetastableStates[X] correspond to the list of relative rotations that correspond to a
         * metastable state/region X. Each rotation is represented by a quaternion.
         * @symmetryQuaternions list of rotation that represents the structural symmetries of the patchy particle if
         * needed to define the states.
         * @param spherePart pointer to equal area spherical partition class with relevant functions
         * @param quaternionPart pointer to volumetric spherical partition class to discretize quaternion space.
         * @param startIndexTransitionStates index of first transition states, if 10, only 9 bound states are supported.
         * If 100, then 99 bound states are suported and so on. This parameter has to be consistent with the one used
         * by the msmrd integrator (see msmrdintegrator for example).
         * @param angularStates number of angular states for discretization
         * @ param tolerance is the maximum distance away from metastable state to still be considered metastable.
         * @param prevsample keeps calue of previous sample when sampling discrete trajectory, useful for
         * coreMSM approach
         */

        patchyDimer(unsigned long Nparticles, int bufferSize);

        void sampleDiscreteTrajectory(double time, std::vector<particle> &particleList) override;

        void setMetastableRegions();

        int getBoundState(quaternion<double> q1, quaternion<double> q2);

    };

}
