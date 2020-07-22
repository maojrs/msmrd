//
// Created by maojrs on 9/13/18.
//

#pragma once
#include "boundary.hpp"


namespace msmrd {
    /**
     * Class for box boundary
     */
    class box : public boundary {
    private:
        double xx = 0.0;
        double yy = 0.0;
        double zz = 0.0;
        std::array<vec3<double>, 3> const normals = {{vec3<double>{-1, 0, 0},
                                                      vec3<double>{0, -1, 0},
                                                      vec3<double>{0, 0, -1}}};

        void enforcePeriodicBoundary(particle &part) override;

        void enforceReflectiveBoundary(particle &part) override;

        void enforceOpenBoundary(particle &part) override;

        void enforcePeriodicBoundary(particleCompound &part) override;

        void enforceReflectiveBoundary(particleCompound &part) override;

        void enforceOpenBoundary(particleCompound &part) override;

    public:
        /**
         * @param xx box length in x-axis from -x/2 to x/2
         * @param yy box length in y-axis from -y/2 to y2
         * @param zz box length in z-axis from -z/2 to z/2
         * @param boundarytype can be periodic, reflective or open
         * @param normals array of the three normals of the planes that constitute the box
         * @param boxsize vector of (xx,yy,zz)
         */
        box(double xxnew, double yynew, double zznew, std::string boundarytype);

        box(vec3<double> newboxsize, std::string boundarytype);

    };

}
