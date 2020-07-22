//
// Created by maojrs on 9/13/18.
//

#pragma once
#include "boundary.hpp"


namespace msmrd{
    /**
     * Class for spherical boundary
     */
    class sphere : public boundary {
    private:
        void enforcePeriodicBoundary(particle &part) override;

        void enforceReflectiveBoundary(particle &part) override;

        void enforceOpenBoundary(particle &part) override;

        void enforcePeriodicBoundary(particleCompound &part) override;

        void enforceReflectiveBoundary(particleCompound &part) override;

        void enforceOpenBoundary(particleCompound &part) override;

    public:
        /**
         * @param radius radius of spherical boundary
         */
        sphere(double radius, std::string boundarytype);

        double getRadius() const { return radius; }

    };

}
