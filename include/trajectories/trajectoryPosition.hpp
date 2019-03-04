//
// Created by maojrs on 3/4/19.
//

#pragma once
#include "trajectory.hpp"

namespace msmrd{
    /**
     * Class to store position only trajectories
     */
    class trajectoryPosition: public trajectory {
    public:

        trajectoryPosition(unsigned long Nparticles, int bufferSize);

        void sample(double time, std::vector<particle> &particleList) override;

        void sampleRelative(double time, std::vector<particle> &particleList) override;

        // Empty function to be overwritten by child classes if neccesary
        void sampleDiscreteTrajectory(double time, std::vector<particle> &particleList) override {};

    };

}