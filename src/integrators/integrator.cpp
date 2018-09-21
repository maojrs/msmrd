//
// Created by dibakma on 27.06.18.
//
#include <array>
#include <utility>
#include "integrators/integrator.hpp"
#include "particle.hpp"

namespace msmrd {
    /**
     * Implementation of integrator abstract class inherited by all child classes
     * @param dt time step
     * @param seed variable for random number generation (Note seed = -1 corresponds to random device)
     * @param randg random number generator based in mt19937
     */
    integrator::integrator(double dt, long seed, std::string bodytype, bool rotation)
            : dt(dt), seed(seed), bodytype(bodytype), rotation(rotation) {
        randg.setSeed(seed);
        clock = 0;
        if (bodytype != "point" && bodytype != "rod" && bodytype != "rigidbody") {
            throw std::runtime_error("Unknown particle bodytype; it should be either point, rod or rigidbody.");
        }
     };

    // Integrate list of particles (need to override in case of MS particles)
    void integrator::integrate(std::vector<particle> &parts) {
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


    // Calculates force and torque from an external potential for point, rod-like and rigidsolid particles
    std::array<vec3<double>, 2> integrator::getExternalForceTorque(particle &part) {
        if (externalPotentialActive) {
            if (bodytype == "point") {
                return externalPot->forceTorque(part.position);
            } else if (bodytype == "rod") {
                return externalRodPot->forceTorque(part.position, part.orientvector);
            } else if (bodytype == "rigidbody") {
                return externalRigidBodyPot->forceTorque(part.position, part.orientation);
            } else {
                throw std::runtime_error("Unknown particle bodytype; it should be either point, rod or rigidbody.");
            };
        }
            // Return zero values if external potential has not been yet set
        else {
            return {vec3<double>(0, 0, 0), vec3<double>(0, 0, 0)};
        }
    };

    // Calculates force and torque due to pair interactions for point, rod-like and rigidsolid particles
    std::array<vec3<double>, 2> integrator::getPairsForceTorque(int partIndex, std::vector<particle> &parts) {
        if (pairPotentialActive) {
            std::array<vec3<double>, 2> forctorq;
            vec3<double> force = vec3<double>(0., 0., 0.);
            vec3<double> torque = vec3<double>(0., 0., 0.);
            if (bodytype == "point") {
                for (int i = 0; i < parts.size(); i++) {
                    if (i != partIndex) {
                        forctorq = pairPot->forceTorque(parts[partIndex].position, parts[i].position);
                        force += forctorq[0];
                        torque += forctorq[1];
                    }
                }
                return {force, torque};
            } else if (bodytype == "rod") {
                for (int i = 0; i < parts.size(); i++) {
                    if (i != partIndex) {
                        forctorq = pairRodPot->forceTorque(parts[partIndex].position, parts[i].position,
                                                           parts[partIndex].orientvector, parts[i].orientvector);
                        force += forctorq[0];
                        torque += forctorq[1];
                    }
                }
                return {force, torque};
            } else if (bodytype == "rigidbody") {
                        for (int i=0; i<parts.size(); i++) {
                            if (i != partIndex) {

                                forctorq = pairRigidBodyPot->forceTorque(parts[partIndex].position, parts[i].position,
                                                                         parts[partIndex].orientation, parts[i].orientation);
                                force += forctorq[0];
                                torque += forctorq[1];
                            }
                        }
                        return {force, torque};
            } else {
                throw std::runtime_error("Unknown particle bodytype. it should be either point, rod or rigidbody.");
            };
        }
            // Return zero values if pair potential has not been yet set
        else {
            return {vec3<double>(0, 0, 0), vec3<double>(0, 0, 0)};
        }
    };


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




