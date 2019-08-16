//
// Created by maojrs on 1/22/19.
//

#pragma once
#include <memory>
#include "trajectories/trajectoryPositionOrientation.hpp"
#include "discretizations/positionOrientationPartition.hpp"
#include "tools.hpp"
#include "H5Cpp.h"

namespace msmrd {
    /**
     * Trajectory class for patchy  dimer trajectories. Patchy dimer here refers to two
     * patchy particles with two patches. The main difference with the general class
     * trajectoryPositionOrientation is that this class can sample the discrete trajectory that is
     * specific for the patchy dimer application. In short words,  it chooses how to discretize the full
     * trajectory of two particle into a discretized trajectory to be analyzed and extracted into
     * a Markov state model. It also implements functionality to discretize trajectories directly loaded
     * from a python array.
     *
     * Patchy dimer will have one unbound state (0), two (or more) bound states (1,2) and
     * angularStates = numSections(numSections+1)/2 angular discrete states (3-3+angularStates).
     * Note in the current indexing implementation only up to 9 bound states are supported.
     */
    class patchyDimer : public trajectoryPositionOrientation {
    private:
        std::unique_ptr<positionOrientationPartition> positionOrientationPart;
        std::array< std::tuple<vec3<double>, quaternion<double>>, 8> boundStates{};
        int maxNumberBoundStates = 10;
        double rLowerBound = 1.25; //1.3; //1.4; //1.25; # 1.4 and 2.2 was a good combination;
        double rUpperBound = 2.25;  //2.5; //2.2  //2.2;   # 1.25 and 2.2 the original configuration. 1.4 and 2.4 the latest
        double tolerancePosition = 0.15;
        double toleranceOrientation = 0.15*2*M_PI;
        int prevsample = 0;
    public:
        /*
         * @positionOrientationPart full six dimensional partition of phase space of relative position and orientation.
         * This is required to sample the discrete trajectory in the transition regions given by this discretization.
         * IMPORTANT: The MSMRD integrator must use the same partition.
         * @boundStates vector of tuples. Each tuple contains a vector and a quaternion indicating one of the 8 bound
         * states. The vector corresponds to the relative position (pos2-pos1) in the frame of reference of particle 1
         * (smaller index) between particle 1 and 2. The fixed frame of reference also assumes particle 1 is in its
         * default initial orientation. The quaternion corresponds to the relative orientation between the two
         * particles, also measured from the fixed frame of reference of particle 1. Thesw bound states are
         * calculated by the setMetastableRegions function.
         * @param maxNumberBoundStates maximum number of bound states supported. It is used to determine how to
         * count (index) the transition states. The state maxNumberBoundStates + 1 will correspond not to a bound state
         * but to the first transition state. This parameter has to be consistent with the one used
         * by the msmrd integrator and the msmrdMarkovModel.
         * @param rLowerBound any relative distance smaller than this value will no longer assign states using
         * the positionOrientationPartition. Instead in the region bounded by r<rLowerBound, either a bound state
         * is assigned or the Core MSM approach will determine the state opf the discrete trajectory, i.e. the previous
         * state will be sampled in the discrete trajectory until a new bound state is reached or r>=rLowerBound.
         * Note rLowerBound < rUpperBound (positionOrientationPart->relativeDistanceCutOff).
         * @param rUpperBound is the upper relative distance limit for the positionOrientationPartition
         * (positionOrientationPart->relativeDistanceCutOff). Any relative distance above will yield the unbound
         * state 0, instead of a transition state in the region rLowerBound <= r < rUpperBound.
         * IMPORTANT: The MSMRD integrator must have the this same value as the relativeDistanceCutOff.
         * @param tolerancePosition is the maximum acceptable difference between the relative position and the
         * calculated relative position of a metstable region to still be considered part of a bound state.
         * @param toleranceOrientation is the maximum acceptable angle-distance difference between the relative
         * orientation and the relative orientation calculated of a given metastable region to still be considerer
         * part of a bound state.
         * @param prevsample keeps calue of previous sample when sampling discrete trajectory, useful for
         * CoreMSM approach. The CoreMSM approach chooses how to discretize the region r<rLowerBound that
         * is not a bound state. CoreMSM uses the value of the previous known bound or transition state until a new
         * bound or transition state is reached.
         */

        patchyDimer(unsigned long Nparticles, int bufferSize);

        patchyDimer(unsigned long Nparticles, int bufferSize, double rLowerBound, double rUpperBound);

        void sampleDiscreteTrajectory(double time, std::vector<particle> &particleList) override;

        int sampleDiscreteState(particle part1, particle part2);

        int getBoundState(vec3<double> relativePosition, quaternion<double> relativeOrientation);

        void setMetastableRegions();



        // Next fucntions are mostly only used when interacting with python or by pybind.

        std::vector<double> discretizeTrajectory(std::vector<std::vector<double>> trajectory);

        // Load H5 directly and discretizes it
        std::vector<double> discretizeTrajectoryH5(std::string filename);

        int getState(particle part1, particle part2);


    };

}
