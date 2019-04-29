//
// Created by maojrs on 2/6/19.
//

#include "integrators/msmrdIntegrator.hpp"


namespace msmrd {

    /**
     * Constructor for MSM/RD integration class, can take discrete time MSM (msm) or a continous-time MSM.
     */
    template<>
    msmrdIntegrator<ctmsm>::msmrdIntegrator(double dt, long seed, std::string particlesbodytype,
                                            double relativeDistanceCutOff, std::vector<ctmsm> MSMlist,
                                            msmrdMSM markovModel, fullPartition positionOrientationPart) :
            overdampedLangevinMarkovSwitch(MSMlist, dt, seed, particlesbodytype), markovModel(markovModel),
            relativeDistanceCutOff(relativeDistanceCutOff), positionOrientationPart(positionOrientationPart) {};

    template<>
    msmrdIntegrator<msm>::msmrdIntegrator(double dt, long seed, std::string particlesbodytype,
                                          double relativeDistanceCutOff, std::vector<msm> MSMlist,
                                          msmrdMSM markovModel, fullPartition positionOrientationPart) :
            overdampedLangevinMarkovSwitch(MSMlist, dt, seed, particlesbodytype), markovModel(markovModel),
            relativeDistanceCutOff(relativeDistanceCutOff), positionOrientationPart(positionOrientationPart) { };



    // Computes possible transitions to bound state and saves them in the event manager. Used by integrate function.
    template<>
    void msmrdIntegrator<ctmsm>::computeTransitions2BoundStates(std::vector<particleMS> &parts) {
        vec3<double> relativePosition;
        quaternion<double> relativeOrientation;
        quaternion<double> refQuaternion;
        int currentTransitionState;
        double transitionTime;
        int nextState;
        for (int i = 0; i < parts.size(); i++) {
            for (int j = i + 1; j < parts.size(); j++) {
                relativePosition = msmrdtools::calculateRelativePosition(parts[i].nextPosition, parts[j].nextPosition,
                                                                         boundaryActive, domainBoundary->getBoundaryType(),
                                                                         domainBoundary->boxsize);

                if (relativePosition.norm() <= relativeDistanceCutOff) {
                    relativeOrientation = parts[j].nextOrientation * parts[i].nextOrientation.conj();
                    refQuaternion = parts[i].nextOrientation.conj();
                    currentTransitionState = positionOrientationPart.getSectionNumber(relativePosition,
                                                                                      relativeOrientation,
                                                                                      refQuaternion);
                    auto transition = markovModel.computeTransition2BoundState(currentTransitionState);
                    transitionTime = std::get<0>(transition);
                    nextState = std::get<1>(transition);
                    eventMgr.addEvent(transitionTime, nextState, i, j, "in");
                }
            }
        }
    }

    /* Computes possible transitions from bound states to to transition states and saves them in the event manager.
     * Used by integrate function. */
    template<>
    void msmrdIntegrator<ctmsm>::computeTransitions2UnboundStates(std::vector<particleMS> &parts) {
        double transitionTime;
        int nextState;
        for (int i = 0; i < parts.size(); i++) {
            // Check if particle[i] is bound (boundTo >= 0) and if bound pairs are only counted once (boundTo > i)
            if (parts[i].boundTo > i) {
                auto transition = markovModel.computeTransition2UnboundState(parts[i].state);
                transitionTime = std::get<0>(transition);
                nextState = std::get<1>(transition);
                eventMgr.addEvent(transitionTime, nextState, i, parts[i].boundTo, "in");
            }
        }
    }



    /* Main integrate function */
    template<>
    void msmrdIntegrator<ctmsm>::integrate(std::vector<particleMS> &parts) {


        /* Integrate only active particles and save next positions/orientations in parts[i].next***.
         * Non-active particles will usually correspond to bound particles */
        for (int i = 0; i < parts.size(); i++) {
            if (parts[i].isActive()) {
                if (parts[i].activeMSM) {
                    integrateOneMS(i, parts, dt);
                } else {
                    integrateOne(i, parts, dt);
                }
            }
        }

        // Enforce boundary and set new positions into parts[i].nextPosition, only if particle is active.
        for (auto &part : parts) {
            if (part.isActive()) {
                if (boundaryActive) {
                    domainBoundary->enforceBoundary(part);
                }
            }
        }

        // Check for transitions to bound states from particles sufficiently close to each other.
        computeTransitions2BoundStates(parts);

        // Check for transitions to unbound states from particles bound to each other.
        computeTransitions2UnboundStates(parts);

        // STILL MISSING A BUNCH OF THINGS

//                    parts[i].deactivate();
//                    parts[j].deactivate();
//                    particleMS boundParticle = particleMS();
//                    boundParticles.push_back();

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