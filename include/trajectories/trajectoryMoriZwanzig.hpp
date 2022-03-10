//
// Created by maojrs on 6/14/21.
//

#pragma once
#include "trajectoryPosition.hpp"
#include "trajectoryPositionVelocity.hpp"

namespace msmrd {
    /**
     * Class to store position only trajectories and auxiliary variable. Trajectory for implementation
     * of the Mori-Zwanzig stochastic closure.
     */
    class trajectoryMoriZwanzig : public trajectoryPositionDistinguished {
    public:
        using trajectoryPositionDistinguished::trajectoryPositionDistinguished;

        void sample(double time, std::vector<particle> &particleList) override;

    };

    /*
     * Same as trajectoryMoriZwanzig but including velocity
     */
    class trajectoryMoriZwanzigVelocity : public trajectoryPositionVelocityDistinguished {
    public:
        using trajectoryPositionVelocityDistinguished::trajectoryPositionVelocityDistinguished;

        void sample(double time, std::vector<particle> &particleList) override;

    };

    /*
     * Same as trajectoryMoriZwanzigVelocity but including an additional aux variable (two aux variables total)
     */
    class trajectoryMoriZwanzigVelocity2 : public trajectoryPositionVelocityDistinguished {
    public:
        using trajectoryPositionVelocityDistinguished::trajectoryPositionVelocityDistinguished;

        void sample(double time, std::vector<particle> &particleList) override;

    };
}
