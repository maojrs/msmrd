//
// Created by maojrs on 2/6/19.
//

#include "integrators/msmrdIntegrator.hpp"


namespace msmrd {

    /**
     * Implementation of base class for MSM/RD integration... to be filled
     */
    msmrdIntegrator::msmrdIntegrator(double dt, long seed, std::string particlesbodytype) :
            overdampedLangevin(dt, seed, particlesbodytype) {

    }


    // Integrate one particle from the particle list main routine (visible only inside the class)
    void msmrdIntegrator::integrateOne(int partIndex, std::vector<particle> &parts, double timestep) {
        vec3<double> force = {0.0, 0.0, 0.0};
        vec3<double> torque = {0.0, 0.0, 0.0};
        translate(parts[partIndex], force, timestep);
        if (rotation) {
            rotate(parts[partIndex], torque, timestep);
        }
    }

    // Main integrate function
    void msmrdIntegrator::integrate(std::vector<particle> &parts, msmrdMSM &masterMSM) {


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

        // Apply MSM/RD coupling for particles sufficiently close to each other.
        vec3<double> relativePosition;
        for (int i = 0; i < parts.size(); i++) {
            for (int j = i + 1; j < parts.size(); j++) {
                relativePosition = parts[j].nextPosition - parts[i].nextPosition;
                if (relativePosition.norm() >= 1.1 && relativePosition.norm() <= 2.1) {
                    auto a = 1; // WAITING FOR IMPLEMENTATION OF DISCRETE MSM TO CHOOSE BETWEEN THE TWO
                }
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



}