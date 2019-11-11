//
// Created by dibakma on 27.06.18.
//
#include <array>
#include <utility>
#include "integrators/integrator.hpp"
#include "particle.hpp"
#include "potentials/potentials.hpp"


namespace msmrd {
    //
    /**
     * Implementation of integrator abstract class inherited by all child classes, note some of its methods
     * are templates and therefore are implemented in the header file
     * @param dt time step
     * @param seed variable for random number generation (Note seed = -1 corresponds to random device).
     * @param particlesbodytype body type of particles to integrate. It determines rotation integrator behavior, can
     * be either point, rod or rigidbody, and it is determined by orientational degrees of freedom, points
     * have no orientation, rods need only one vector and rigidsolid requires a complete quaternion).
     * @param randg random number generator based in mt19937
     */
    integrator::integrator(double dt, long seed, std::string particlesbodytype)
            : dt(dt), seed(seed), particlesbodytype(particlesbodytype) {
        randg.setSeed(seed);
        clock = 0;
        forceField.resize(0);
        torqueField.resize(0);
        if (particlesbodytype != "point") {
            rotation = true;
        }
        if (particlesbodytype != "point" && particlesbodytype != "rod" && particlesbodytype != "rigidbody") {
            throw std::invalid_argument("Unknown particles bodytype; it should be either point, rod, rigidbody");
        }
    };

    // Integrate list of particles (need to override in case of MS particles)
    void integrator::integrate(std::vector<particle> &parts) {
        vec3<double> force;
        vec3<double> torque;

        // Calculate forces and torques and save them into forceField and torqueField
        calculateForceTorqueFields(parts);
        // Integrate and save next positions/orientations in parts[i].next***
        for (int i = 0; i < parts.size(); i++) {
            integrateOne(i, parts, dt);
        }

        // Enforce boundary; sets new positions into parts[i].nextPosition (only if particle is active)
        enforceBoundary(parts);

        // Update position/orientation based on parts[i].nextPosition, parts[i].nextOrientation
        updatePositionOrientation(parts);

        // Updates time
        clock += dt;
    }

    // Calculates relative distance (p2-p1) of two vectors, p1, p2, taking into account possible periodic boundary
    vec3<double> integrator::calculateRelativePosition(vec3<double> p1, vec3<double> p2) {
        // Calculate relative distance. If box periodic boundary, take that into account.
        if (boundaryActive and domainBoundary->getBoundaryType() == "periodic") {
            auto boxsize = domainBoundary->boxsize;
            return msmrdtools::distancePeriodicBox(p1, p2, boxsize);
        } else {
            return p2 - p1;
        }
    }


    // Incorporates custom boundary into integrator
    void integrator::setBoundary(boundary *bndry) {
        boundaryActive = true;
        domainBoundary = bndry;
        if (pairPotentialActive) {
            pairPot->setBoundary(bndry);
        }
    }


    /* Integrator set potential functions for external potentials (for particles of type particle or its children
     * like particleMS). */
    void integrator::setExternalPotential(externalPotential *pot) {
        externalPotentialActive = true;
        externalPot = pot;
    }

    /* Integrator set potential functions for pair potentials (for particles of type particle or its children
     * like particleMS). */
    void integrator::setPairPotential(pairPotential *pot) {
        pairPotentialActive = true;
        pairPot = pot;
        if (boundaryActive) {
            pairPot->setBoundary( domainBoundary );
        }
    }

}




