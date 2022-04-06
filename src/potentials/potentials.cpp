//
// Created by maojrs on 8/21/18.
//
#include "potentials/potentials.hpp"

namespace msmrd {
    /**
     *  Implementation of non-abstract functions of externalPotential abstract class (default constructor in header)
     *  These functions are needed for PyBinding forceTorque functions, since pyBind only vector output.
     */


    /* Force Torque for binding functions */
    std::vector<std::vector<double>> externalPotential::forceTorquePyBind(particle &part) {
        std::array<vec3<double>, 2> forceTorquex = forceTorque(part);
        return msmrdtools::array2Dtovec2D(forceTorquex);
    }

    std::vector<std::vector<double>> pairPotential::forceTorquePyBind(particle &part1, particle &part2) {
        std::array<vec3<double>, 4> forceTorquex = forceTorque(part1, part2);
        return msmrdtools::array2Dtovec2D(forceTorquex);
    }

    // Incorporates integrator's boundary into potential
    void pairPotential::setBoundary(boundary *bndry) {
        boundaryActive = true;
        domainBoundary = bndry;
    }

    // Calculates relative distance (p2-p1) of two vectors, p1, p2, taking into account possible periodic boundary
    vec3<double> pairPotential::relativePosition(const vec3<double> p1, const vec3<double> p2) {
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


    /*
     * Implementations for combined potentials
     */

    void combinedPairPotential::addPotential(std::shared_ptr<pairPotential> pairPotentialPtr) {
        potentials.push_back(std::move(pairPotentialPtr));
    }

    double combinedPairPotential::evaluate(particle &part1, particle &part2) {
        double potentialVal = 0.0;
        for(const auto& potential: potentials) {
            potentialVal += potential->evaluate(part1, part2);
        }
        return potentialVal;
    };

    std::array<vec3<double>, 4> combinedPairPotential::forceTorque(particle &part1, particle &part2) {
        vec3<double> force = vec3<double>(0, 0, 0);
        vec3<double> torque = vec3<double>(0, 0, 0);
        for(const auto& potential: potentials) {
            auto forctorq = potential->forceTorque(part1, part2);
            force += forctorq[0];
            torque += forctorq[1];
        }
        return {force, torque, -1.0*force, -1.0*torque};
    };


    void combinedExternalPotential::addPotential(std::shared_ptr<externalPotential> externalPotentialPtr) {
        potentials.push_back(std::move(externalPotentialPtr));
    }

    double combinedExternalPotential::evaluate(particle &part1) {
        double potentialVal = 0.0;
        for(const auto& potential: potentials) {
            potentialVal += potential->evaluate(part1);
        }
        return potentialVal;
    };

    std::array<vec3<double>, 2> combinedExternalPotential::forceTorque(particle &part1) {
        vec3<double> force = vec3<double>(0, 0, 0);
        vec3<double> torque = vec3<double>(0, 0, 0);
        for(const auto& potential: potentials) {
            auto forctorq = potential->forceTorque(part1);
            force += forctorq[0];
            torque += forctorq[1];
        }
        return {force, torque};
    };

}