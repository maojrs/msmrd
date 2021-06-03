//
// Created by maojrs on 3/4/19.
//

#pragma once
#include "trajectory.hpp"

namespace msmrd{
    /**
     * Class to store trajectories with position and orientation (given by a quaternion)
     */
    class trajectoryPositionOrientation : public trajectory {
    public:

        trajectoryPositionOrientation(unsigned long Nparticles, int bufferSize);

        void sample(double time, std::vector<particle> &particleList) override;

        void sampleRelative(double time, std::vector<particle> &particleList) override;

        // Empty function to be overwritten by child classes if neccesary
        void sampleDiscreteTrajectory(double time, std::vector<particle> &particleList) override {};

        void printTime();
    };

    class trajectoryPositionOrientationType : public trajectoryPositionOrientation {
    public:
        using trajectoryPositionOrientation::trajectoryPositionOrientation;

        void sample(double time, std::vector<particle> &particleList) override;
    };


    /**
     * Class to store trajectories with position, orientation and state
     */
    class trajectoryPositionOrientationState : public trajectoryPositionOrientation {
    public:
        using trajectoryPositionOrientation::trajectoryPositionOrientation;

        void sample(double time, std::vector<particle> &particleList) override;

    };


    class trajectoryPositionOrientationStateType : public trajectoryPositionOrientationState {
    public:
        using trajectoryPositionOrientationState::trajectoryPositionOrientationState;

        void sample(double time, std::vector<particle> &particleList) override;
    };

}
