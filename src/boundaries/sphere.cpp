//
// Created by maojrs on 9/13/18.
//

#include "boundaries/sphere.hpp"


namespace msmrd {
    /**
     * Box boundary class implementation
     * @param x box length in x-axis from -x/2 to x/2
     * @param y box length in y-axis from -y/2 to y2
     * @param z box length in z-axis from -z/2 to z/2
     */
    sphere::sphere(double newradius, std::string boundarytype) : boundary(boundarytype) {
        radius = newradius;
    };

    // Enforces periodic box boundary condition
    void sphere::enforcePeriodicBoundary(particle &part) {
        throw std::runtime_error("Periodic boundary condition not available for spherical boundary.");
    };

    // Enforces reflective box boundary condition
    void sphere::enforceReflectiveBoundary(particle &part) {
        if (part.nextPosition.norm() > radius && part.position.norm() <=radius ) {
            vec3<double> r0 = part.position;
            vec3<double> dr = part.nextPosition - part.position;
            // Find intersection point r0 + al*dr with sphere
            double A = dr * dr;
            double B = 2.0 * (r0 * dr);
            double C = (r0 * r0) - radius * radius;
            double discriminant = B * B - 4.0 * A * C;
            if (discriminant < 0) {
                throw std::range_error("Discriminant is zero, error in spherical reflective boundary. "
                                       "Check all your particles are initially contained in boundary.");
            }
            double al = (-B + std::sqrt(discriminant)) / (2.0 * A); //positive root is the correct one
            vec3<double> intersection = r0 + al * dr;
            // Normal to sphere at intersection point
            vec3<double> normal = -1.0*intersection / intersection.norm();
            // Obtain reflected vector from boundary function and assign into newPosition
            vec3<double> reflectedvec = reflectVector(r0, dr, intersection, normal);
            part.setNextPosition(reflectedvec);
        }
    };

    // Enforces open boundary condition
    void sphere::enforceOpenBoundary(particle &part) {
        if (part.nextPosition.norm() > radius)  {
            part.deactivate();
        }
    }

}