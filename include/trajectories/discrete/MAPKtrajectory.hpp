//
// Created by maojrs on 5/3/21.
//

#pragma once
#include <memory>
#include "trajectories/discrete/discreteTrajectory.hpp"
#include "discretizations/positionOrientvectorPartition.hpp"
#include "tools.hpp"

namespace msmrd {
/**
 * Trajectory class for patchy protein trajectory for the MAPK model (specific application). This class
 * samples the discrete trajectory. We just need to specify how to discretize the full trajectory by
 * setting all the possible bound states between all the particles.
 *
 * Note discretization used here should match the discretization of MSM/RD for consistent results.
 *
 * Also note we are overriding virtual and nonvirtual parts of the parent class discreteTrajectory.
 * Functionality easy to brake if not done carefully.
 */

    /* Patchy protein MAPK trajectory implementation with one unbound state (0)
     * and 4 bound states (1,2,3,4). The bound states 1 and 2 correspond to the kinase binding to sites
     * one and two, while bound states 3 and 4 correspond to the phosphatase binding to the same
     * two binding sites. */
    class MAPKtrajectory : public discreteTrajectory<4> {
    protected:
        std::unique_ptr<positionOrientvectorPartition> positionOrientvectorPart;
        std::array< std::tuple<vec3<double>, vec3<double>>, 4> boundStates{};
        /*
         * @positionOrientvectorPart five dimensional partition of phase space of relative position and
         * orientation given by an orientvector. This is required to sample the discrete trajectory in
         * the transition regions given by this discretization. IMPORTANT: The MSMRD integrator must use
         * the same partition. Also this pointer should be set in constructor. This substitutes
         * positionOrientationPart functionality of the parent discreteTrajectory class.
         * @boundStates vector of tuples. Each tuple contains two vectors indicating each one of the
         * bound states. The vector corresponds to the relative position (pos2-pos1) in the frame of reference of
         * particle 1 (main particle). The fixed frame of reference also assumes particle 1
         * is in its default initial orientation. The second vector corresponds to the orientation vector
         * of the second particle in the frame of reference of particle 1. These bound states are
         * calculated by the setBoundStates function. This substitutes the variable with the same
         * name in the parent discreteTrajectory class.
         */

        public:

        MAPKtrajectory(unsigned long Nparticles, int bufferSize);

        MAPKtrajectory(unsigned long Nparticles, int bufferSize, double rLowerBound, double rUpperBound);

        void setBoundStates();

        int sampleDiscreteState(particle part1, particle part2) override;

        int getBoundState(vec3<double> relativePosition, vec3<double> orientVector);

        vec3<double> getRelativePosition(int boundStateIndex);

        quaternion<double> getRelativeOrientation(int boundStateIndex);

        };


}