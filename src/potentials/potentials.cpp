//
// Created by maojrs on 8/21/18.
//
#include <tuple>
#include "potentials/potentials.hpp"
#include "tools.hpp"
#include "vec3.hpp"

namespace msmrd {
    /**
     *  Implementation of non-abstract functions of externalPotential abstract class (default constructor in header)
     *  These functions are needed for PyBinding forceTorque functions, since pyBind only vector output.
     */


    /* Force Torque for binding functions */
    std::vector<std::vector<double>> externalPotential::forceTorquePyBind(const particle &part) {
        std::array<vec3<double>, 2> forceTorquex = forceTorque(part);
        return msmrdtools::array2Dtovec2D(forceTorquex);
    }

    std::vector<std::vector<double>> pairPotential::forceTorquePyBind(const particle &part1, const particle &part2) {
        std::array<vec3<double>, 4> forceTorquex = forceTorque(part1, part2);
        return msmrdtools::array2Dtovec2D(forceTorquex);
    }

    // Incorporates integrator's boundary into potential
    void pairPotential::setBoundary(boundary *bndry) {
        boundaryActive = true;
        domainBoundary = bndry;
    }

    // Calculates relative distance (p2-p1) of two vectors, p1, p2, taking into account possible periodic boundary
    vec3<double> pairPotential::relativePosition(vec3<double> p1, vec3<double> p2) {
        // Calculate relative distance. If box periodic boundary, take that into account.
        if (boundaryActive and domainBoundary->getBoundaryType() == "periodic") {
            auto boxsize = domainBoundary->boxsize;
            return msmrdtools::distancePeriodicBox(p1, p2, boxsize);
        } else {
            return p2 - p1;
        }
    };

    /* Calculates relative distance (p2-p1) of two vectors, p1, p2, taking into account possible periodic boundary and
     * returns virtual position of p1 as well as the relative distance */
    std::array<vec3<double>, 2> pairPotential::relativePositionComplete(vec3<double> p1, vec3<double> p2){
        // Calculate relative distance. If box periodic boundary, take that into account.
        if (boundaryActive and domainBoundary->getBoundaryType() == "periodic") {
            vec3<double> boxsize = domainBoundary->boxsize;
            return msmrdtools::distancePeriodicBoxComplete(p1, p2, boxsize);
        } else {
            return {p1, p2 - p1};
        }
    };

}