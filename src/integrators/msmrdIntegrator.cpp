//
// Created by maojrs on 2/6/19.
//

#include "integrators/msmrdIntegrator.hpp"


namespace msmrd {

    // Template constructors in header


    /* Computes possible transitions to bound states from particles sufficiently close to each other (in transition
     * states) and saves them in the event manager. Used by integrate function. */
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
                // Need special function to calculate relative position, in case we use a periodic boundary.
                if (boundaryActive) {
                    relativePosition = msmrdtools::calculateRelativePosition(parts[i].nextPosition, parts[j].nextPosition,
                            boundaryActive, domainBoundary->getBoundaryType(), domainBoundary->getBoxsize());
                } else {
                    relativePosition = parts[j].nextPosition - parts[i].nextPosition;
                }

                if (relativePosition.norm() < relativeDistanceCutOff) {
                    relativeOrientation = parts[j].nextOrientation * parts[i].nextOrientation.conj();
                    refQuaternion = parts[i].nextOrientation.conj();
                    currentTransitionState = positionOrientationPart->getSectionNumber(relativePosition,
                                                                                      relativeOrientation,
                                                                                      refQuaternion);
                    // PROBLEM ON LINE ABOVE OR BELOW ...
                    auto transition = markovModel.computeTransition2BoundState(currentTransitionState);
                    transitionTime = std::get<0>(transition);
                    nextState = std::get<1>(transition);
                    eventMgr.addEvent(transitionTime, nextState, i, j, "in");
                }
            }
        }
    }

    /* Computes possible transitions from unbound states to transition states and saves them in the event manager.
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

    /* Makes particles with indexes iIndex and jIndex in the particle list transition to a bound state. Note
     * always iIndex < jIndex should hold. Also particle with smaller index is the ones that remains active
     * to model position/orientation of bound complex. */
    template<>
    void msmrdIntegrator<ctmsm>::transition2BoundState(std::vector<particleMS> &parts, int iIndex,
                                                       int jIndex, int endState) {
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
        // Average bound particle position and orientation (save on particle with smaller index).
        parts[iIndex].nextPosition = 0.5*(parts[iIndex].position + parts[jIndex].position);
        parts[iIndex].nextOrientation = msmrdtools::quaternionSlerp(parts[iIndex].orientation,
                                                                parts[jIndex].orientation, 0.5);
        // Assign constant distant position to inactive particle ( could be useful for visualization).
        parts[jIndex].nextPosition = {10000000.0, 10000000.0, 10000000.0};
    }


    /* Makes particles with indexes iIndex and jIndex in the particle list transition to an unbound state. Note
     * always iIndex < jIndex should hold.*/
    template<>
    void msmrdIntegrator<ctmsm>::transition2UnboundState(std::vector<particleMS> &parts, int iIndex,
                                                       int jIndex, int endState) {
        int iNewState;
        int jNewState;
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
        auto iDiff = MSMlist[iPartType].Dlist[iNewState];
        auto iDiffRot = MSMlist[iPartType].Drotlist[iNewState];
        auto jDiff = MSMlist[jPartType].Dlist[jNewState];
        auto jDiffRot = MSMlist[jPartType].Drotlist[jNewState];
        parts[iIndex].setDs(iDiff, iDiffRot);
        parts[jIndex].setDs(jDiff, jDiffRot);
        // Extract section intervals in partition corresponding to the endState.
        auto sections = positionOrientationPart->getSectionIntervals(endState);
        auto phiInterval = std::get<0>(sections); //polar
        auto thetaInterval = std::get<1>(sections); //azimuthal
        auto quatRadInterval = std::get<2>(sections);
        auto quatPhiInterval = std::get<3>(sections); //polar
        auto quatThetaInterval = std::get<4>(sections); //azimuthal
        // Calculate new relative positions and orientations I AM HERE NOW, STILL NEED TO IMPLEMENT THIS
        auto relPosition = relativeDistanceCutOff * randg.uniformSphereSection(phiInterval, thetaInterval);
        // Calculate realtive orientation from discretization (partition)
        auto randomQuat = randg.uniformShellSection(quatRadInterval, quatPhiInterval, quatThetaInterval);
        double sQuat = std::sqrt(1 - randomQuat.norm());
        // Recover quaternion from its representation as a vector inside the unit 3D sphere
        quaternion<double> relOrientation = {sQuat, randomQuat};
        // Set next positions and orientations based on the relative ones (parts[iIndex] keeps track of position)
        parts[iIndex].nextPosition = parts[iIndex].position - 0.5*relPosition;
        parts[jIndex].nextPosition = parts[iIndex].nextPosition + relPosition;
        // More complicate approach for orientation is possible but likely uneccesary.
        parts[iIndex].nextOrientation = 1.0 * parts[iIndex].orientation;
        parts[jIndex].nextOrientation = relOrientation * parts[iIndex].nextOrientation;
    }

    /* Removes unrealized events where the particles drifted a distance apart beyond the relativeDistanceCutOff.
     * This can be more optimally included inside the computeTransitions2BoundStates function. However, the
     * code is more clear if left separate. Better left for future code optimization. */
    template<>
    void msmrdIntegrator<ctmsm>::removeUnrealizedEvents(std::vector<particleMS> &parts) {
        int numEvents = eventMgr.getNumEvents();
        for (int i = numEvents - 1; i>=0; --i) {
            auto event = eventMgr.getEvent(i);
            auto iIndex = std::get<2>(event)[0];
            auto jIndex = std::get<2>(event)[1];
            vec3<double> relativePosition;
            if (boundaryActive) {
                relativePosition = msmrdtools::calculateRelativePosition(parts[iIndex].nextPosition,
                        parts[jIndex].nextPosition, boundaryActive, domainBoundary->getBoundaryType(),
                                                                         domainBoundary->getBoxsize());
            } else{
                relativePosition = parts[jIndex].nextPosition - parts[iIndex].nextPosition;
            }
            // Remove event if particles drifted apart
            if (relativePosition.norm() >= relativeDistanceCutOff) {
                eventMgr.removeEvent(i);
            }
        }
    }





    /* Main integrate function */
    template<>
    void msmrdIntegrator<ctmsm>::integrate(std::vector<particleMS> &parts) {

        // Calculate forces and torques and save them into forceField and torqueField (needed even if there are zero)
        calculateForceTorqueFields<particleMS>(parts); // Should only run once if forces are zero.

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

        /* Remove unrealized previous events. Also compute transitions to bound states and to unbound states and add
         * them as new events in the event manager. Finally resort event list by descending waiting time (last in
         * list happen first) */
        removeUnrealizedEvents(parts);
        computeTransitions2BoundStates(parts); // PROBLEM WITH THIS FUNCTION!!! CHECK
        computeTransitions2UnboundStates(parts);
        eventMgr.sortDescending();


        // Check for events that should happen in this time step and make them happen.
        int numEvents = eventMgr.getNumEvents();
        for (int i = numEvents - 1; i>=0; --i) {
            // Break loop if only future events are left (eventTime>0)
            if (eventMgr.getEventTime(i) > 0) {
                break;
            } else {
                // Load event data
                auto event = eventMgr.getEvent(i);
                auto residualTime = std::get<0>(event);
                auto endState = std::get<1>(event);
                auto iIndex = std::get<2>(event)[0];
                auto jIndex = std::get<2>(event)[1];
                auto inORout = std::get<3>(event);
                // Make event happen
                if (inORout == "in") {
                    transition2BoundState(parts, iIndex, jIndex, endState);
                } else if (inORout == "out") {
                    transition2UnboundState(parts, iIndex, jIndex, endState);
                }
                // Remove event from event list
                eventMgr.removeEvent(i);
            }
        }

        // Enforce boundary and set new positions into parts[i].nextPosition (only if particle is active).
        for (auto &part : parts) {
            if (part.isActive() and boundaryActive) {
                domainBoundary->enforceBoundary(part);
            }
        }

        /* Update positions and orientations (sets calculated next position/orientation
         * calculated by integrator and boundary as current position/orientation). Note states
         * are modified directly and don't need to be updated. */
        for (int i = 0; i < parts.size(); i++) {
            if (parts[i].isActive()) {
                parts[i].updatePosition();
                if (rotation) {
                    parts[i].updateOrientation();
                }
            }
        }

        // Advance global time and time in event manager
        clock += dt;
        eventMgr.advanceTime(dt);
    }


}