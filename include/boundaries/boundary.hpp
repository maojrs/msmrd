//
// Created by maojrs on 9/4/18.
//
#pragma once
#include <iostream>
#include "particle.hpp"
#include "vec3.hpp"


namespace msmrd {
    /**
     * Abstract base class to set the domain boundary
     */
    class boundary {
    protected:
        std::string boundarytype;
       /*
        * @param boundarytype can be periodic, reflective or open
        */

        virtual void enforcePeriodicBoundary(particle &part) = 0;

        virtual void enforceReflectiveBoundary(particle &part) = 0;

        virtual void enforceOpenBoundary(particle &part) = 0;

    public:
        boundary(std::string boundarytype);

        void enforceBoundary(particle &part);
    };

}
