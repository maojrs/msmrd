//
// Created by maojrs on 6/14/21.
//

#pragma once
#include "trajectoryPosition.hpp"

namespace msmrd {
    /**
     * Class to store position only trajectories
     */
    class trajectoryMoriZwanzig : public trajectoryPositionDistinguished {
    public:
        using trajectoryPositionDistinguished::trajectoryPositionDistinguished;

        void sample(double time, std::vector<particle> &particleList) override;

    };
}
