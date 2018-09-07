//
// Created by maojrs on 9/5/18.
//
#include <array>
#include "boundary.hpp"
#include "particle.hpp"


namespace msmrd {
    /**
     * All boundary classes implementations
     */


    /**
     * Box boundary class implementation
     * @param x box length in x-axis from -x/2 to x/2
     * @param y box length in y-axis from -y/2 to y2
     * @param z box length in z-axis from -z/2 to z/2
     */
    box::box(double xx, double yy, double zz, std::string boundarytype) : xx(xx), yy(yy), zz(zz),
                                                                          boundarytype(boundarytype) {
        boxsize = vec3<double>{xx, yy, zz};
    };

    box::box(vec3<double> boxsize, std::string boundarytype) : boxsize(boxsize), boundarytype(boundarytype) {
        xx = boxsize[0];
        yy = boxsize[1];
        zz = boxsize[2];
    };


    /* Main function to enforce boundary, note enforce boundary acts on "nextPositions" of the particles, so
     * it can keep better track of the previous position when the boundary requires so. */
    void box::enforceBoundary(particle &part) {
        if (boundarytype == "periodic") {
            enforcePeriodicBoundary(part);
        } else if (boundarytype == "reflective") {
            enforceReflectiveBoundary(part);
        } else if (boundarytype == "open") {
            enforceOpenBoundary(part);
        } else {
            throw std::runtime_error("Unknown boundary type; it should be either periodic, reflective or open.");
        }
    };

    // Enforces periodic box boundary condition
    void box::enforcePeriodicBoundary(particle &part) {
        for (int i = 0; i < 3; i++) {
            if (part.nextPosition[i] >= boxsize[i] / 2) { part.nextPosition[i] -= boxsize[i]; }
            if (part.nextPosition[i] <= -boxsize[i] / 2) { part.nextPosition[i] += boxsize[i]; }
        }
    };

    // Enforces reflective box boundary condition
    void box::enforceReflectiveBoundary(particle &part) {
        vec3<double> reflectedVec;
        vec3<double> invnormal;
        double scalarProd;
    // TO BE IMPLEMENTED
    //    for (int i = 0; i<3; i++) {
    //        if (part.nextPosition[i] >= boxsize[i] / 2) {
    //            scalarProd = normals[i]*part.nextPosition;
    //            reflectedVec = part.nextPosition - 2*normals[i]*scalarProd;
    //        }
    //        if (part.nextPosition[i] <= -boxsize[i] / 2) {
    //            invnormal = -1*normals[i];
    //            scalarProd = invnormal*part.nextPosition;
    //            reflectedVec = part.nextPosition - 2*invnormal*scalarProd;
    //        }
    //    }
    };

    void box::enforceOpenBoundary(particle &part) {
        if (std::abs(part.nextPosition[0]) > xx ||
            std::abs(part.nextPosition[1]) > yy ||
            std::abs(part.nextPosition[2]) > zz) {
            part.deactivate();
        }
    }

}