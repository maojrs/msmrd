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
     * @param seed variable for random number generation (Note seed = -1 corresponds to random device)
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
            throw std::runtime_error("Unknown particles bodytype; it should be either point, rod, rigidbody");
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

        // Enforce boundary and set new positions into parts[i].nextPosition
        for (auto &part : parts) {
            if (boundaryActive) {
                domainBoundary->enforceBoundary(part);
            }
        }
        /* Update positions and orientations (sets calculated next position/orientation
         * calculated by integrator and boundary as current position/orientation). */
        for (auto &part : parts) {
            part.updatePosition();
            if (rotation) {
                part.updateOrientation();
            }
        }
        clock += dt;
    }

    // Incorporates custom boundary into integrator
    void integrator::setBoundary(boundary *bndry) {
        boundaryActive = true;
        domainBoundary = bndry;
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
    }

}




