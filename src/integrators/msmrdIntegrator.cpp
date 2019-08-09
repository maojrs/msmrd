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
                /* Only compute new transition if particles drifted into transition region for
                 * the first time, i.e. empty event */
                auto previousEvent = eventMgr.getEvent(i, j);
                if (previousEvent.inORout == "empty") {
                    // Need special function to calculate relative position, in case we use a periodic boundary.
                    relativePosition = calculateRelativePosition(parts[i].nextPosition, parts[j].nextPosition);
                    if (relativePosition.norm() < relativeDistanceCutOff) {
                        if (rotation) {
                            relativeOrientation = parts[i].nextOrientation.conj() * parts[j].nextOrientation;
                            refQuaternion = parts[i].nextOrientation.conj();
                            currentTransitionState = positionOrientationPart->getSectionNumber(relativePosition,
                                                                                               relativeOrientation,
                                                                                               refQuaternion);
                        } else {
                            currentTransitionState = positionPart->getSectionNumber(relativePosition);
                        }
                        auto transition = markovModel.computeTransition2BoundState(currentTransitionState);
                        transitionTime = std::get<0>(transition);
                        nextState = std::get<1>(transition);
                        eventMgr.addEvent(transitionTime, i, j, currentTransitionState, nextState, "in");

                    }
                }
            }
        }
    }

    /* Computes possible transitions from bound states to transition states and saves them in the event manager.
     * Used by integrate function. */
    template<>
    void msmrdIntegrator<ctmsm>::computeTransitions2UnboundStates(std::vector<particleMS> &parts) {
        double transitionTime;
        int nextState;
        std::tuple<double, int> transition;
        for (int i = 0; i < parts.size(); i++) {
            // Only compute transition if particle is bound to another particle.
            // If particle[i] is bound (boundTo > 0); if bound pairs are only counted once (boundTo > i)
            if (parts[i].boundTo > i) {
                /* Only compute transition if particles drifted into bound state for
                 * the first time, i.e. empty event */
                auto previousEvent = eventMgr.getEvent(i, parts[i].boundTo);
                if (previousEvent.inORout == "empty") {
                    transition = markovModel.computeTransition2UnboundState(parts[i].state);
                    transitionTime = std::get<0>(transition);
                    nextState = std::get<1>(transition);
                    eventMgr.addEvent(transitionTime, i, parts[i].boundTo, parts[i].state, nextState, "out");
                }
            }
        }
    }

    /* Computes possible transitions from bound states to other bound states and saves them in the event manager.
     * Used by integrate function. */
    template<>
    void msmrdIntegrator<ctmsm>::computeTransitionsBetweenBoundStates(std::vector<particleMS> &parts) {
        double transitionTime;
        int nextState;
        std::tuple<double, int> transition;
        for (int i = 0; i < parts.size(); i++) {
            // Check if particle[i] is bound (boundTo > 0); if bound pairs are only counted once (boundTo > i)
            if (parts[i].boundTo > i) {
                auto previousEvent = eventMgr.getEvent(i, parts[i].boundTo);
                /* Only compute transition if particles drifted into new bound state for
                 * the first time, i.e. empty event */
                if (previousEvent.inORout == "empty") {
                    transition = markovModel.computeTransitionBetweenBoundStates(parts[i].state);
                    transitionTime = std::get<0>(transition);
                    nextState = std::get<1>(transition);
                    eventMgr.addEvent(transitionTime, i, parts[i].boundTo, parts[i].state, nextState, "inside");
                }
            }
        }
    }

    /* Makes particles with indexes iIndex and jIndex in the particle list transition to a bound state. Note
     * always iIndex < jIndex should hold. Also particle with smaller index is the one that remains active
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
        // Set diffusion coefficients in bound state (note states start counting from 1, not zero)
        parts[iIndex].setDs(markovModel.Dboundlist[endState-1], markovModel.DboundRotlist[endState-1]);
        // Average bound particle position and orientation (save on particle with smaller index).
        if (rotation) {
            parts[iIndex].nextOrientation = msmrdtools::quaternionSlerp(parts[iIndex].orientation,
                                                                        parts[jIndex].orientation, 0.5);
        }
        parts[iIndex].nextPosition = 0.5*(parts[iIndex].position + parts[jIndex].position);
        // Assign constant distant position to inactive particle ( could be useful for visualization).
        parts[jIndex].nextPosition = {10000000.0, 10000000.0, 10000000.0};
        parts[jIndex].updatePosition();
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
        std::array<double, 2> phiInterval;
        std::array<double, 2> thetaInterval;
        if (not rotation) {
            auto sections = positionPart->getAngles(endState);
            phiInterval = std::get<0>(sections); //polar
            thetaInterval = std::get<1>(sections); //azimuthal
        } else {
            std::array<double, 2> quatRadInterval;
            std::array<double, 2> quatPhiInterval;
            std::array<double, 2> quatThetaInterval;
            auto sections = positionOrientationPart->getSectionIntervals(endState);
            phiInterval = std::get<0>(sections); //polar
            thetaInterval = std::get<1>(sections); //azimuthal
            quatRadInterval = std::get<2>(sections);
            quatPhiInterval = std::get<3>(sections); //polar
            quatThetaInterval = std::get<4>(sections); //azimuthal
            // Calculate relative orientation from discretization (partition)
            auto randomQuat = randg.uniformShellSection(quatRadInterval, quatPhiInterval, quatThetaInterval);
            double sQuat = std::sqrt(1 - randomQuat.norm());
            // Recover quaternion from its representation as a vector inside the unit 3D sphere
            quaternion<double> relOrientation = {sQuat, randomQuat};
            // More complicate approach for orientation is possible but likely uneccesary.
            parts[iIndex].nextOrientation = 1.0 * parts[iIndex].orientation;
            parts[jIndex].nextOrientation = relOrientation * parts[iIndex].nextOrientation;
        }

        // Calculate new relative positions and orientations
        auto relPosition = relativeDistanceCutOff * randg.uniformSphereSection(phiInterval, thetaInterval);

        // Set next positions and orientations based on the relative ones (parts[iIndex] keeps track of position)
        parts[iIndex].nextPosition = parts[iIndex].position - 0.5*relPosition;
        parts[jIndex].nextPosition = parts[iIndex].nextPosition + relPosition;
    }


    /* Transitions already bound particles with indexes iIndex and jIndex in the particle list
     * to another bound state. */
    template<>
    void msmrdIntegrator<ctmsm>::transitionBetweenBoundStates(std::vector<particleMS> &parts, int iIndex,
                                                       int jIndex, int endState) {
        // Set state for particle
        parts[iIndex].setState(endState);
        parts[jIndex].setState(endState);
        // Set diffusion coefficients in bound state (note states start counting from 1, not zero)
        parts[iIndex].setDs(markovModel.Dboundlist[endState-1], markovModel.DboundRotlist[endState-1]);
    }


    /* Removes unrealized events where unbound particles drifted a distance apart beyond the relativeDistanceCutOff,
     * or when zero rates yielded infinite values. This can be more optimally included inside the
     * computeTransitions2BoundStates function. However, the code is more clear if left separate. Better left for
     * future code optimization. */
    template<>
    void msmrdIntegrator<ctmsm>::removeUnrealizedEvents(std::vector<particleMS> &parts) {
        // Important not to use range loop, iterator loop better since events are being erased.
        auto it = eventMgr.eventDictionary.begin();
        while(it != eventMgr.eventDictionary.end()) {
            auto transitionTime = it->second.waitTime;
            auto iIndex = it->second.part1Index;
            auto jIndex = it->second.part2Index;
            vec3<double> relativePosition;
            // Remove event if transition time is infinity
            if (std::isinf(transitionTime)) {
                //eventMgr.removeEvent(iIndex, jIndex);
                //erase() will return the next iterator
                it = eventMgr.eventDictionary.erase(it);
                continue;
            }
            // If particles in unbound state and relative position is larger than cutOff, remove event.
            if (parts[iIndex].boundTo == -1 and parts[jIndex].boundTo == -1) {

                relativePosition = calculateRelativePosition(parts[iIndex].nextPosition, parts[jIndex].nextPosition);

                // Remove event if particles drifted apart
                if (relativePosition.norm() >= relativeDistanceCutOff) {
                    //eventMgr.removeEvent(iIndex, jIndex);
                    it = eventMgr.eventDictionary.erase(it);
                    continue;
                }
            }
            it++;
        }
    }





    /* Main integrate function */
    template<>
    void msmrdIntegrator<ctmsm>::integrate(std::vector<particleMS> &parts) {

        /* Calculate forces and torques and save them into forceField and torqueField. For the MSM/RD this will
         * in general zero, so only needs to be run once. */
        if (firstrun) {
            calculateForceTorqueFields<particleMS>(parts);
            firstrun = false;
        }

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
        computeTransitions2BoundStates(parts);
        computeTransitions2UnboundStates(parts);
        computeTransitionsBetweenBoundStates(parts);


        // Advance global time and in event manager (make events in [t,t+dt) happen)
        clock += dt;
        eventMgr.advanceTime(dt);

        // Check for events that should happen in this time step and make them happen.
        for (auto &thisEvent : eventMgr.eventDictionary) {
            auto transitionTime = thisEvent.second.waitTime;

            // Break loop if only future events are left (eventTime>0)
            if (transitionTime < 0) {
                // Load event data
                auto endState = thisEvent.second.endState;
                auto iIndex = thisEvent.second.part1Index;
                auto jIndex = thisEvent.second.part2Index;
                auto inORout = thisEvent.second.inORout;
                // Make event happen
                if (thisEvent.second.inORout == "in") {
                    transition2BoundState(parts, iIndex, jIndex, endState);
                } else if (inORout == "out") {
                    transition2UnboundState(parts, iIndex, jIndex, endState);
                } else if (inORout == "inside") {
                    transitionBetweenBoundStates(parts, iIndex, jIndex, endState);
                }
                // Remove event from event list
                eventMgr.removeEvent(iIndex, jIndex);
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

    }


}