//
// Created by maojrs on 9/4/18.
//
#pragma once
#include <iostream>
#include "particle.hpp"
#include "particleCompound.hpp"
#include "vec3.hpp"


namespace msmrd {
    /**
     * Abstract base class to set the domain boundary
     */
    class boundary {
    protected:
        std::string boundarytype;
       /*
        * @param boundarytype can be periodic, reflective, open or none
        * @param radius size of boundary in case of spherical boundary (simplifies code to have it in parent class)
        * @param boxsize boundary size in case of box boundary (simplifies code to have it in parent class)
        */

        virtual void enforcePeriodicBoundary(particle &part) = 0;

        virtual void enforceReflectiveBoundary(particle &part) = 0;

        virtual void enforceOpenBoundary(particle &part) = 0;

        virtual void enforcePeriodicBoundary(particleCompound &part) = 0;

        virtual void enforceReflectiveBoundary(particleCompound &part) = 0;

        virtual void enforceOpenBoundary(particleCompound &part) = 0;

    public:

        double radius = 0;

        vec3<double> boxsize = {0, 0, 0};

        boundary(std::string boundarytype);

        template< typename PARTICLE >
        void enforceBoundary(PARTICLE &particle);

        std::string getBoundaryType() const { return boundarytype; }

        vec3<double> getBoxsize() const { return boxsize; }


        /** Function to reflect vector
         * @param dr vector pointing from previous position to current position
         * @param r0 vector of previous position
         * @param intersection cooridnates of intersection with boundary
         * @param normal unitary inner normal in boundary at intersection point
         * returns corresponding portion of vector dr that was reflected in boundary
         */
         vec3<double> reflectVector(vec3<double> r0, vec3<double> dr, vec3<double> intersection, vec3<double> normal);
    };


    /* Main template function to enforce boundary, note enforce boundary acts on "nextPositions" of the particles, so
     * it can keep better track of the previous position when the boundary requires so. */
    template <typename PARTICLE>
    void boundary::enforceBoundary(PARTICLE &part) {
        if (boundarytype == "periodic") {
            enforcePeriodicBoundary(part);
        } else if (boundarytype == "reflective") {
            enforceReflectiveBoundary(part);
        } else if (boundarytype == "open") {
            enforceOpenBoundary(part);
        } else if (boundarytype != "none") {
            throw std::runtime_error("Unknown boundary type; it should be either periodic, reflective or open.");
        }
    }

}
