//
// Created by maojrs on 5/29/19.
//

#pragma once
#include "boundary.hpp"


namespace msmrd {
    /**
     * Class for inactive boundary (sometimes necessary for simplified functionality of msmrdIntegrator).
     * This is the default class used by integrator when boundary is neither set nor active.
     */
    class noBoundary : public boundary {
    private:
        // Empty functions since there is no actual boundary

        void enforcePeriodicBoundary(particle &part) override {};

        void enforceReflectiveBoundary(particle &part) override {};

        void enforceOpenBoundary(particle &part) override {};

        void enforcePeriodicBoundary(particleCompound &part) override {};

        void enforceReflectiveBoundary(particleCompound &part) override {};

        void enforceOpenBoundary(particleCompound &part) override {};

    public:

        noBoundary() : boundary("none") {};
    };

}
