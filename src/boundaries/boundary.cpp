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

}