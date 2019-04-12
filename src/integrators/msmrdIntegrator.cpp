//
// Created by maojrs on 2/6/19.
//

#include "integrators/msmrdIntegrator.hpp"


namespace msmrd {

    /**
     * Implementation of base class for MSM/RD integration... to be filled
     */
    msmrdIntegrator::msmrdIntegrator(double dt, long seed, std::string particlesbodytype, msmrdMSM markovModel,
                                     fullPartition positionOrientationPart) :
            overdampedLangevin(dt, seed, particlesbodytype), markovModel(markovModel),
            positionOrientationPart(positionOrientationPart) { };


    /* Integrate one particle from the particle list main routine (visible only inside the class), only used
     * if particles are far away. */
    void msmrdIntegrator::integrateOne(int partIndex, std::vector<particle> &parts, double timestep) {
        vec3<double> force = {0.0, 0.0, 0.0};
        vec3<double> torque = {0.0, 0.0, 0.0};
        if ( parts[partIndex].isActive() ) {
            translate(parts[partIndex], force, timestep);
            if (rotation) {
                rotate(parts[partIndex], torque, timestep);
            }
        }
    }

    // Main integrate function
    void msmrdIntegrator::integrate(std::vector<particle> &parts, msmrdMSM &masterMSM) {


        // Integrate only active particles and save next positions/orientations in parts[i].next***
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
        quaternion<double> relativeOrientation;
        quaternion<double> refQuaternion;
        int currentTransitionState;
        double transitionTime;
        int nextState;
        for (int i = 0; i < parts.size(); i++) {
            for (int j = i + 1; j < parts.size(); j++) {
                relativePosition = calculateRelativePosition(parts[i], parts[j]);

                if (relativePosition.norm() <= 2.2) {
                    relativeOrientation = parts[j].nextOrientation * parts[i].nextOrientation.conj();
                    refQuaternion = parts[i].nextOrientation.conj();
                    currentTransitionState = positionOrientationPart.getSectionNumber(relativePosition,
                                                                                      relativeOrientation,
                                                                                      refQuaternion);
                    auto transition = markovModel.computeTransition(currentTransitionState);
                    transitionTime = std::get<0>(transition);
                    nextState = std::get<1>(transition);
                    // STILL MISSING A BUNCH OF THINGS
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


    vec3<double> msmrdIntegrator::calculateRelativePosition(particle &p1, particle &p2){
        vec3<double> relPosition;
        if (boundaryActive and domainBoundary->getBoundaryType() == "periodic") {
            auto boxsize = domainBoundary->boxsize;
            relPosition = msmrdtools::distancePeriodicBox(p2.nextPosition, p1.nextPosition, boxsize);
        } else {
            relPosition = p2.nextPosition - p1.nextPosition;
        }
        return relPosition;
    }

}