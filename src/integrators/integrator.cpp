//
// Created by dibakma on 27.06.18.
//
#include <array>
#include <utility>
#include "integrators/integrator.hpp"
#include "particle.hpp"

namespace msmrd {
    //
    /**
     * Implementation of integrator abstract class inherited by all child classes, note some of its methods
     * are templates and therefore are implemented in the header file
     * @param dt time step
     * @param seed variable for random number generation (Note seed = -1 corresponds to random device)
     * @param randg random number generator based in mt19937
     */
    integrator::integrator(double dt, long seed, std::string particlestype, bool rotation)
            : dt(dt), seed(seed), particlestype(particlestype), rotation(rotation) {
        randg.setSeed(seed);
        clock = 0;
        forceField.resize(0);
        torqueField.resize(0);
        if (particlestype != "point" && particlestype != "rod" && particlestype != "rigidbody" &&
            particlestype != "pointMS" && particlestype != "rodMS" && particlestype != "rigidbodyMS") {
            throw std::runtime_error("Unknown particle bodytype; it should be either point, rod, rigidbody,"
                                     "pointMS, rodMS or rigibodyMS.");
        }
     };

    // Integrate list of particles (need to override in case of MS particles)
    void integrator::integrate(std::vector<particle> &parts) {
        vec3<double> force;
        vec3<double> torque;

        // Calculate forces and torques and save them into forceField and torqueField
        calculateForceTorqueFields<particle>(parts);
        // Integrate and save next positions/orientations in parts[i].next***
        for (int i = 0; i < parts.size(); i++) {
            integrateOne(i, parts, dt);
        }

        // Enforce boundary and set new positions into parts[i].nextPosition
        for (int i = 0; i < parts.size(); i++) {
            if (boundaryActive) {
                domainBoundary->enforceBoundary(parts[i]);
            }
        }
        /* Update positions and orientations (sets calculated next position/orientation
         * calculated by integrator and boundary as current position/orientation). */
        for (int i = 0; i < parts.size(); i++) {
                parts[i].updatePosition();
            if (rotation) {
                parts[i].updateOrientation();
            }
        }
        clock += dt;
    }

    // Incorporates custom boundary into integrator
    void integrator::setBoundary(boundary *bndry) {
        boundaryActive = true;
        domainBoundary = bndry;
    }

    // Incorporates custom external potential functions into integrator
    void integrator::setExternalPotential(externalPotential<> *pot) {
        externalPotentialActive = true;
        externalPot = pot;
    }

    // Incorporates custom external potential functions for rod-like particles into integrator
    void integrator::setExternalRodPotential(externalPotential<vec3<double>> *pot) {
        externalPotentialActive = true;
        externalRodPot = pot;
    }

    // Incorporates custom external potential functions for rigidbody particles into integrator
    void integrator::setExternalRigidBodyPotential(externalPotential<quaternion<double>> *pot) {
        externalPotentialActive = true;
        externalRigidBodyPot = pot;
    }

    // Incorporates custom pair potential functions into integrator
    void integrator::setPairPotential(pairPotential<> *pot) {
        pairPotentialActive = true;
        pairPot = pot;
    }

    // Incorporates custom pair potential function for rod-like particles into integrator
    void integrator::setPairRodPotential(pairPotential<vec3<double>, vec3<double>> *pot) {
        pairPotentialActive = true;
        pairRodPot = pot;
    }

    // Incorporates custom pair potential function for rigidbody particles into integrator
    void integrator::setPairRigidBodyPotential(pairPotential<quaternion<double>, quaternion<double>> *pot) {
        pairPotentialActive = true;
        pairRigidBodyPot = pot;
    }

}




