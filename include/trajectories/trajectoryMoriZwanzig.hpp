//
// Created by maojrs on 6/14/21.
//

#pragma once
#include "trajectoryPosition.hpp"

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
}
