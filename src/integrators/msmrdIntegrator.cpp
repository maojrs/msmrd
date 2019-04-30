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
                                            int numParticleTypes, double relativeDistanceCutOff,
                                            std::vector<ctmsm> MSMlist, msmrdMSM markovModel,
                                            fullPartition positionOrientationPart) :
            overdampedLangevinMarkovSwitch(MSMlist, dt, seed, particlesbodytype), markovModel(markovModel),
            numParticleTypes(numParticleTypes), relativeDistanceCutOff(relativeDistanceCutOff),
            positionOrientationPart(positionOrientationPart) {
        if (MSMlist.size() != numParticleTypes) {
            std::__throw_range_error("Number of MSMs provided in MSMlist should match number of unbound particle"
                                     "types in the simulation");
        }
    };

    template<>
    msmrdIntegrator<msm>::msmrdIntegrator(double dt, long seed, std::string particlesbodytype,
                                          int numParticleTypes,  double relativeDistanceCutOff,
                                          std::vector<msm> MSMlist, msmrdMSM markovModel,
                                          fullPartition positionOrientationPart) :
            overdampedLangevinMarkovSwitch(MSMlist, dt, seed, particlesbodytype), markovModel(markovModel),
            numParticleTypes(numParticleTypes), relativeDistanceCutOff(relativeDistanceCutOff),
            positionOrientationPart(positionOrientationPart) {
        if (MSMlist.size() != numParticleTypes) {
            std::__throw_range_error("Number of MSMs provided in MSMlist should match number of unbound particle"
                                     "types in the simulation");
        }
    };



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
                eventMgr.addEvent(transitionTime, nextState, i, parts[i].boundTo, "out");
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

        // Sort list, check for events that should happen and make them happen.
        eventMgr.sort();
        int numEvents = eventMgr.getNumEvents();
        int eventCounter = 0;
        int iNewState;
        int jNewState;
        for (int i = 0; i < numEvents; i++) {
            if (eventMgr.getEventTime(i) > 0) {
                break;
            }
            else {
                eventCounter += 1;
                auto event = eventMgr.getEvent(i);
                auto residualTime = std::get<0>(event);
                auto endState = std::get<1>(event);
                auto iIndex = std::get<2>(event)[0];
                auto jIndex = std::get<2>(event)[1];
                auto inORout = std::get<3>(event);
                // Make transition to bound state happen (smaller index remains the active particle)
                if (inORout == "in") {
                    /* Establish pair connection in particle class and deactivate particle with larger index ( only
                     * particle with smaller index remains active to represent the movement of the bound particle) */
                    parts[iIndex].boundTo = jIndex;
                    parts[jIndex].boundTo = iIndex;
                    parts[jIndex].deactivate();
                    // Set state for particle
                    parts[iIndex].setState(endState);
                    parts[jIndex].setState(endState);
                    // Deactivate unbound MSM behavior if applicable
                    parts[iIndex].activeMSM = false;
                    parts[jIndex].activeMSM = false;
                    // Set diffusion coefficients in bound state
                    parts[iIndex].setDs(markovModel.Dboundlist[endState], markovModel.DboundRotlist[endState]);
                    /* Average bound particle position and orientation (save on particle with smaller index),
                     * send deactivated particle far away ( could be useful for visualization)*/
                    parts[iIndex].position = 0.5*(parts[iIndex].position + parts[jIndex].position);
                    parts[iIndex].orientation = msmrdtools::quaternionSlerp(parts[iIndex].orientation,
                                                                            parts[jIndex].orientation, 0.5);
                    parts[jIndex].position = {10000000.0, 10000000.0, 10000000.0};

                }
                //Make transition to unbound state happen
                else if (inORout == "out") {
                    // Eliminate pair connection by resetting to default value -1 and activate particles.
                    parts[iIndex].boundTo = -1;
                    parts[jIndex].boundTo = -1;
                    parts[iIndex].activate();
                    parts[jIndex].activate();
                    /* Calculates next states and activates MSM if there is a transition matrix (size larger than one).
                     * More complex behavior to calculate new states possible */
                    iNewState = 0;
                    jNewState = 0;
                    auto iPartType = parts[iIndex].type;
                    auto jPartType = parts[jIndex].type;
                    if (MSMlist[iPartType].tmatrix.size() > 1) {
                        parts[iIndex].activeMSM = true;
                        iNewState = randg.uniformInteger(0, MSMlist[iPartType].tmatrix.size());
                    }
                    if (MSMlist[jPartType].tmatrix.size() > 1) {
                        parts[jIndex].activeMSM = true;
                        jNewState = randg.uniformInteger(0, MSMlist[jPartType].tmatrix.size());
                    }
                    // Set new states
                    parts[iIndex].setState(iNewState);
                    parts[jIndex].setState(jNewState);
                    // Set diffusion coefficients
                    parts[iIndex].setDs(MSMlist[iPartType].Dlist[iNewState],
                                        MSMlist[iPartType].Drotlist[iNewState]);
                    parts[jIndex].setDs(MSMlist[jPartType].Dlist[jNewState],
                                        MSMlist[jPartType].Drotlist[jNewState]);
                    // Sets new positions and orientations I AM HERE NOW, STILL NEED TO IMPLEMENT THIS
                }
            }
        }

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

        // Advance global time and time in event manager
        clock += dt;
        eventMgr.advanceTime(dt);
    }





}