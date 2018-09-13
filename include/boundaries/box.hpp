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
        double xx, yy, zz;
        std::array<vec3<double>, 3> normals = {{vec3<double>{-1, 0, 0},
                                                vec3<double>{0, -1, 0},
                                                vec3<double>{0, 0, -1}}};

        void enforcePeriodicBoundary(particle &part) override;

        void enforceReflectiveBoundary(particle &part) override;

        void enforceOpenBoundary(particle &part) override;

    public:
        vec3<double> boxsize;

        /**
         * @param xx box length in x-axis from -x/2 to x/2
         * @param yy box length in y-axis from -y/2 to y2
         * @param zz box length in z-axis from -z/2 to z/2
         * @param boundarytype can be periodic, reflective or open
         * @param normals array of the three normals of the planes that constitute the box
         * @param boxsize vector of (xx,yy,zz)
         */
        box(double xx, double yy, double zz, std::string boundarytype);

        box(vec3<double> boxsize, std::string boundarytype);

        std::string getBoundaryType() const { return boundarytype; }
    };

}
