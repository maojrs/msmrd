//
// Created by maojrs on 9/13/18.
//

#include "boundaries/box.hpp"


namespace msmrd {
    /**
     * Box boundary class implementation
     * @param x box length in x-axis from -x/2 to x/2
     * @param y box length in y-axis from -y/2 to y/2
     * @param z box length in z-axis from -z/2 to z/2
     */
    box::box(double xx, double yy, double zz, std::string boundarytype)
            : xx(xx), yy(yy), zz(zz), boundary(boundarytype) {
        boxsize = vec3<double>{xx, yy, zz};
    };

    box::box(vec3<double> newboxsize, std::string boundarytype) : boundary(boundarytype) {
        boxsize = newboxsize;
        xx = boxsize[0];
        yy = boxsize[1];
        zz = boxsize[2];
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
        //throw std::runtime_error("Reflective boundary for box is not yet implemented.");
        for (int i = 0; i<3; i++) {
            if (part.nextPosition[i] >= boxsize[i] / 2) {
                vec3<double> r0 = part.nextPosition;
                vec3<double> dr = part.nextPosition - part.position;
                // Find intersection point
                double al = (boxsize[i] / 2 - r0[i])/dr[i];
                vec3<double> intersection = r0 + al*dr;
                // Obtain reflected vector from boundary function and assign into newPosition
                vec3<double> reflectedvec = reflectVector(r0, dr, intersection, normals[i]);
                part.setPosition(intersection); // Need to asign in case it bounces in a corner (other loop iteration)
                part.setNextPosition(reflectedvec);
            }
            if (part.nextPosition[i] <= -boxsize[i] / 2) {
                invnormal = -1*normals[i];
                vec3<double> r0 = part.nextPosition;
                vec3<double> dr = part.nextPosition - part.position;
                // Find intersection point
                double al = (-boxsize[i] / 2 - r0[i])/dr[i];
                vec3<double> intersection = r0 + al*dr;
                // Obtain reflected vector from boundary function and assign into newPosition
                vec3<double> reflectedvec = reflectVector(r0, dr, intersection, invnormal);
                part.setPosition(intersection); // Need to asign in case it bounces in a corner (other loop iteration)
                part.setNextPosition(reflectedvec);
            }
        }
    };

    // Enforces open boundary condition
    void box::enforceOpenBoundary(particle &part) {
        if (std::abs(part.nextPosition[0]) > xx ||
            std::abs(part.nextPosition[1]) > yy ||
            std::abs(part.nextPosition[2]) > zz) {
            part.deactivate();
        }
    }

}