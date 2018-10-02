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
    integrator::integrator(double dt, long seed, std::string particlesbodytype, bool rotation)
            : dt(dt), seed(seed), particlesbodytype(particlesbodytype), rotation(rotation) {
        randg.setSeed(seed);
        clock = 0;
        forceField.resize(0);
        torqueField.resize(0);
        if (particlesbodytype != "point" && particlesbodytype != "rod" && particlesbodytype != "rigidbody" &&
            particlesbodytype != "pointMix" && particlesbodytype != "rodMix" && particlesbodytype != "rigidbodyMix") {
            throw std::runtime_error("Unknown particles bodytype; it should be either point, rod, rigidbody,"
                                     "pointMix, rodMix or rigibodyMix.");
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


    /*
     * Integrator set potential functions for external potentials.
     */


    // Incorporates custom external potential functions into integrator
    void integrator::setExternalPotential(externalPotential<> *pot) {
        if (particlesbodytype != "point") {
            throw std::runtime_error("This potential requires particles bodytype = 'point'. ");
        }
        externalPotentialActive = true;
        externalPot = pot;
    }

    // Incorporates custom external potential functions for mix of point particles into integrator
    void integrator::setExternalMixPotential(externalPotential<int> *pot) {
        if (particlesbodytype != "pointMix") {
            throw std::runtime_error("This potential requires particles bodytype = 'pointMix'. ");
        }
        externalPotentialActive = true;
        externalMixPot = pot;
    }

    // Incorporates custom external potential functions for rod-like particles into integrator
    void integrator::setExternalRodPotential(externalPotential<vec3<double>> *pot) {
        if (particlesbodytype != "rod") {
            throw std::runtime_error("This potential requires particles bodytype = 'rod'. ");
        }
        externalPotentialActive = true;
        externalRodPot = pot;
    }

    // Incorporates custom external potential functions for mix of rod-like particles into integrator
    void integrator::setExternalRodMixPotential(externalPotential<vec3<double>, int> *pot) {
        if (particlesbodytype != "rodMix") {
            throw std::runtime_error("This potential requires particles bodytype = 'rodMix'. ");
        }
        externalPotentialActive = true;
        externalRodMixPot = pot;
    }

    // Incorporates custom external potential functions for rigidbody particles into integrator
    void integrator::setExternalRigidBodyPotential(externalPotential<quaternion<double>> *pot) {
        if (particlesbodytype != "rigidbody") {
            throw std::runtime_error("This potential requires particles bodytype = 'rigidbody'. ");
        }
        externalPotentialActive = true;
        externalRigidBodyPot = pot;
    }

    // Incorporates custom external potential functions for mix of rigidbody particles into integrator
    void integrator::setExternalRigidBodyMixPotential(externalPotential<quaternion<double>, int> *pot) {
        if (particlesbodytype != "rigidbodyMix") {
            throw std::runtime_error("This potential requires particles bodytype = 'rigidbodyMix'. ");
        }
        externalPotentialActive = true;
        externalRigidBodyMixPot = pot;
    }


    /*
     * Integrator set potential functions for pair potentials.
     */


    // Incorporates custom pair potential functions into integrator
    void integrator::setPairPotential(pairPotential<> *pot) {
        if (particlesbodytype != "point") {
            throw std::runtime_error("This potential requires particles bodytype = 'point'. ");
        }
        pairPotentialActive = true;
        pairPot = pot;
    }

    // Incorporates custom pair potential functions for mix of point particles into integrator
    void integrator::setPairMixPotential(pairPotential<int, int> *pot) {
        if (particlesbodytype != "pointMix") {
            throw std::runtime_error("This potential requires particles bodytype = 'pointMix'. ");
        }
        pairPotentialActive = true;
        pairMixPot = pot;
    }

    // Incorporates custom pair potential function for rod-like particles into integrator
    void integrator::setPairRodPotential(pairPotential<vec3<double>, vec3<double>> *pot) {
        if (particlesbodytype != "rod") {
            throw std::runtime_error("This potential requires particles bodytype = 'rod'. ");
        }
        pairPotentialActive = true;
        pairRodPot = pot;
    }

    // Incorporates custom pair potential function for mix of rod-like particles into integrator
    void integrator::setPairRodMixPotential(pairPotential<vec3<double>, vec3<double>, int, int> *pot) {
        if (particlesbodytype != "rodMix") {
            throw std::runtime_error("This potential requires particles bodytype = 'rodMix'. ");
        }
        pairPotentialActive = true;
        pairRodMixPot = pot;
    }

    // Incorporates custom pair potential function for rigidbody particles into integrator
    void integrator::setPairRigidBodyPotential(pairPotential<quaternion<double>, quaternion<double>> *pot) {
        if (particlesbodytype != "rigidbody") {
            throw std::runtime_error("This potential requires particles bodytype = 'rigidbody'. ");
        }
        pairPotentialActive = true;
        pairRigidBodyPot = pot;
    }

    // Incorporates custom pair potential function for mix of rigidbody particles into integrator
    void integrator::setPairRigidBodyMixPotential(pairPotential<quaternion<double>, quaternion<double>, int, int> *pot) {
        if (particlesbodytype != "rigidbodyMix") {
            throw std::runtime_error("This potential requires particles bodytype = 'rigidbodyMix'. ");
        }
        pairPotentialActive = true;
        pairRigidBodyMixPot = pot;
    }

}




