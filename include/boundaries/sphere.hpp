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
        double radius;

        void enforcePeriodicBoundary(particle &part) override;

        void enforceReflectiveBoundary(particle &part) override;

        void enforceOpenBoundary(particle &part) override;

    public:
        /**
         * @param radius radius of spherical boundary
         */
        sphere(double radius, std::string boundarytype);

        std::string getBoundaryType() const { return boundarytype; }
        double getRadius() const { return radius; }

    };

}
