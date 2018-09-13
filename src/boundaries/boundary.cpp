//
// Created by maojrs on 9/5/18.
//
#include "boundaries/boundary.hpp"
#include "particle.hpp"


namespace msmrd {
    /**
     * Common functionality to all boundary classes (boundary parent class) implementation
     */
    boundary::boundary(std::string boundarytype) : boundarytype(boundarytype) {};


    /* Main function to enforce boundary, note enforce boundary acts on "nextPositions" of the particles, so
     * it can keep better track of the previous position when the boundary requires so. */
    void boundary::enforceBoundary(particle &part) {
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

    /** Function to reflect vector
     * @param dr vector pointing from previous position to current position
     * @param r0 vector of previous position
     * @param intersection cooridnates of intersection with boundary
     * @param normal unitary inner normal in boundary at intersection point
     * returns corresponding portion of vector dr that was reflected in boundary
     */
    vec3<double> boundary::reflectVector(vec3<double> r0, vec3<double> dr,
                                         vec3<double> intersection, vec3<double> normal){
        // Calculate reflected vector (reflection of vector in plane formula), and rescale
        vec3<double> reflectedvec = dr - 2.0 * (dr * normal) * normal;
        vec3<double> remaindervec = r0 + dr - intersection;
        reflectedvec = remaindervec.norm()*reflectedvec / reflectedvec.norm();
        // Translate reflected vector
        reflectedvec = intersection + reflectedvec;
        return reflectedvec;
    }

}