//
// Created by maojrs on 8/19/19.
//

#include "integrators/msmrdIntegratorDiscrete.hpp"

namespace msmrd {

    // Template constructors in header

    /* Computes the current transition state in the discretization for an unbound pair of
     * particles. If their relative position is larger than the discretization limit (radialBounds[1]),
     * it returns -1. Otherwise it returns the transition state (in the original discrete trajectories indexing)
     * of the two particles in the discretization.*/
    template<>
    int msmrdIntegratorDiscrete<ctmsm>::computeCurrentTransitionState(particleMS &part1, particleMS &part2) {
        int currentTransitionState = -1;
        vec3<double> relativePosition;
        quaternion<double> relativeOrientation;
        quaternion<double> refQuaternion;
        int index0 = markovModel.getMaxNumberBoundStates();
        // Need special function to calculate relative position, in case we have a periodic boundary.
        relativePosition = calculateRelativePosition(part1.nextPosition, part2.nextPosition);
        if (relativePosition.norm() < radialBounds[1]) {
            if (rotation) {
                relativeOrientation = part1.nextOrientation.conj() * part2.nextOrientation;
                refQuaternion = part1.nextOrientation.conj();
                currentTransitionState = positionOrientationPart->getSectionNumber(relativePosition,
                                                                                   relativeOrientation,
                                                                                   refQuaternion);
            } else {
                currentTransitionState = positionPart->getSectionNumber(relativePosition);
            }
            // Calculate transition with Markov model, note reindexing by index0 required.
            currentTransitionState = index0 + currentTransitionState;
        }
        return currentTransitionState;
    }


    /* Computes possible transitions from transition states to bound states or other transition states from
     * particles sufficiently close to each other (in transition states) and saves them in the event manager.
     * Used by integrate function. */
    template<>
    void msmrdIntegratorDiscrete<ctmsm>::computeTransitionsFromTransitionStates(std::vector<particleMS> &parts) {
        int currentTransitionState;
        double transitionTime;
        int nextState;
        int index0 = markovModel.getMaxNumberBoundStates();
        // Loop over all pairs of particles without repetition with i < j
        for (int i = 0; i < parts.size(); i++) {
            for (int j = i + 1; j < parts.size(); j++) {
                // Only compute transitions if both particles are in unbound state.
                if (parts[i].boundTo == -1 and parts[j].boundTo == -1) {
                    /* Computes new transition if particles drifted into transition region for
                     * the first time, i.e. empty event and relativeDistance < radialBounds[1], or if
                     * particles transitioned between transition states. */
                    auto previousEvent = eventMgr.getEvent(i, j);
                    if (previousEvent.eventType == "empty") { // returns -1 if |relativePosition| > radialBounds[1]
                        currentTransitionState = computeCurrentTransitionState(parts[i], parts[j]);
                    } else if (previousEvent.eventType == "inTransition") {
                        currentTransitionState = previousEvent.endState;
                        eventMgr.removeEvent(i, j);
                    }
                    // If valid currentTransitionState (see computeCurrentTransitionState), calculate next transition.
                    if (currentTransitionState != -1) {
                        auto transition = markovModel.calculateTransition(currentTransitionState);
                        transitionTime = std::get<0>(transition);
                        nextState = std::get<1>(transition);
                        if (nextState <= index0) {
                            eventMgr.addEvent(transitionTime, i, j, currentTransitionState, nextState, "binding");
                        } else {
                            eventMgr.addEvent(transitionTime, i, j, currentTransitionState,
                                              nextState, "transition2transition");
                        }
                    }
                }
            }
        }
    }


    /* Computes possible transitions from bound states to other bound states or unbound states (transition states)
     * and saves them in the event manager. Used by integrate function. */
    template<>
    void msmrdIntegratorDiscrete<ctmsm>::computeTransitionsFromBoundStates(std::vector<particleMS> &parts) {
        double transitionTime;
        int nextState;
        std::tuple<double, int> transition;
        int index0 = markovModel.getMaxNumberBoundStates();
        for (int i = 0; i < parts.size(); i++) {
            // Only compute transition if particle is bound to another particle.
            // If particle[i] is bound, boundTo > 0; if bound pairs are only to be counted once, then boundTo > i
            if (parts[i].boundTo > i) {
                /* Only compute transition if particles switched into a given bound state for
                 * the first time, i.e. empty event */
                auto previousEvent = eventMgr.getEvent(i, parts[i].boundTo);
                if (previousEvent.eventType == "empty") {
                    transition = markovModel.calculateTransition(parts[i].boundState);
                    transitionTime = std::get<0>(transition);
                    nextState = std::get<1>(transition);
                    // Distinguish between events bound to bound transition and unbinding events
                    if (nextState <= index0) {
                        eventMgr.addEvent(transitionTime, i, parts[i].boundTo,
                                          parts[i].state, nextState, "bound2bound");
                    } else {
                        eventMgr.addEvent(transitionTime, i, parts[i].boundTo,
                                          parts[i].state, nextState, "unbinding");
                    }
                }
            }
        }
    }


    /* Makes particles with indexes iIndex and jIndex in the particle list transition to a bound state. Note
     * always iIndex < jIndex should hold. Also particle with smaller index is the one that remains active
     * to model position/orientation of bound complex. */
    template<>
    void msmrdIntegratorDiscrete<ctmsm>::transition2BoundState(std::vector<particleMS> &parts, int iIndex,
                                                       int jIndex, int endState) {
        /* Establish pair connection in particle class and deactivate particle with larger index ( only
         * particle with smaller index remains active to represent the movement of the bound particle) */
        parts[iIndex].boundTo = jIndex;
        parts[jIndex].boundTo = iIndex;
        parts[jIndex].deactivate();
        // Set bound state for particle
        parts[iIndex].setBoundState(endState);
        parts[jIndex].setBoundState(endState);
        // Set diffusion coefficients in bound state (note states start counting from 1, not zero)
        parts[iIndex].setDs(markovModel.Dlist[endState-1], markovModel.Drotlist[endState-1]);
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
    void msmrdIntegratorDiscrete<ctmsm>::transition2UnboundState(std::vector<particleMS> &parts, int iIndex,
                                                         int jIndex, int endState) {

        int iNewState;
        int jNewState;
        // Redefine endstate indexing, so it is understood by the partition/discretization.
        int index0 = markovModel.getMaxNumberBoundStates();
        endState = endState - index0;

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
        /* Set new unbound states (which also eliminates pair connections by resetting boundTo and
         * boundState to -1) and activate particles. */
        parts[iIndex].setState(iNewState);
        parts[jIndex].setState(jNewState);
        parts[iIndex].activate();
        parts[jIndex].activate();
        // Sets diffusion coefficients
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

        // Calculate new relative positions and orientations by sampling uniformly in spherical section.
        auto relPosition = randg.uniformShellSection(radialBounds, phiInterval, thetaInterval);

        // Set next positions and orientations based on the relative ones (parts[iIndex] keeps track of position)
        parts[iIndex].nextPosition = parts[iIndex].position - 0.5*relPosition;
        parts[jIndex].nextPosition = parts[iIndex].nextPosition + relPosition;
    }


    /* Transitions of already bound particles with indexes iIndex and jIndex in the particle list
     * to another bound state. */
    template<>
    void msmrdIntegratorDiscrete<ctmsm>::transitionBetweenBoundStates(std::vector<particleMS> &parts, int iIndex,
                                                              int jIndex, int endState) {
        // Set state for particle
        parts[iIndex].setBoundState(endState);
        parts[jIndex].setBoundState(endState);
        /* Set diffusion coefficients in bound state (although not always neccesary for bound states,
         * we convert to the proper indexing.) */
        int MSMindex = markovModel.getMSMindex(endState);
        parts[iIndex].setDs(markovModel.Dlist[MSMindex], markovModel.Drotlist[MSMindex]);
    }

    /* Transitions of particles in transition region with indexes iIndex and jIndex in the particle list
     * to another transition state. */
    template<>
    void msmrdIntegratorDiscrete<ctmsm>::transitionBetweenTransitionStates(int iIndex, int jIndex) {
        /* Add label to event to indicate to computeTransitionsFromTransitionStates function that
         * a new event needs to be calculated, using previousEvent.endState as the initial state */
        auto previousEvent = eventMgr.getEvent(iIndex, jIndex);
        previousEvent.eventType = "inTransition";
    }


    /* Removes unrealized events where unbound particles drifted a distance apart beyond the upper radial bound,
     * or when zero rates yielded infinite values. This can be more optimally included inside the
     * computeTransitionsFromTransitionStates function. However, the code is more clear if left separate.
     * Better left for future code optimization. */
    template<>
    void msmrdIntegratorDiscrete<ctmsm>::removeUnrealizedEvents(std::vector<particleMS> &parts) {
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
                if (relativePosition.norm() >= radialBounds[1]) {
                    //eventMgr.removeEvent(iIndex, jIndex);
                    it = eventMgr.eventDictionary.erase(it);
                    continue;
                }
            }
            it++;
        }
    }


    /* Apply events in event manager that should happen during the current time step. */
    template<>
    void msmrdIntegratorDiscrete<ctmsm>::applyEvents(std::vector<particleMS> &parts) {
        for (auto &thisEvent : eventMgr.eventDictionary) {
            auto transitionTime = thisEvent.second.waitTime;
            /* Only apply events that should happen in during this
             * timestep (they correspond to transitionTime < 0) */
            if (transitionTime < 0) {
                // Load event data
                auto iIndex = thisEvent.second.part1Index;
                auto jIndex = thisEvent.second.part2Index;
                auto endState = thisEvent.second.endState;
                auto eventType = thisEvent.second.eventType;
                // Make event happen (depending on event type) and remove event once it has happened
                if (eventType == "binding") {
                    transition2BoundState(parts, iIndex, jIndex, endState);
                    eventMgr.removeEvent(iIndex, jIndex);
                } else if (eventType == "unbinding") {
                    transition2UnboundState(parts, iIndex, jIndex, endState);
                    eventMgr.removeEvent(iIndex, jIndex);
                } else if (eventType == "bound2bound") {
                    transitionBetweenBoundStates(parts, iIndex, jIndex, endState);
                    eventMgr.removeEvent(iIndex, jIndex);
                } else if (eventType == "transition2transition") {
                    transitionBetweenTransitionStates(iIndex, jIndex);
                }
            }
        }
    }


    /* Main integrate function */
    template<>
    void msmrdIntegratorDiscrete<ctmsm>::integrate(std::vector<particleMS> &parts) {

        /* Calculate forces and torques and save them into forceField and torqueField. For the MSM/RD this will
         * in general be zero, so only needs to be run once. */
        if (firstrun) {
            calculateForceTorqueFields<particleMS>(parts);
            firstrun = false;
        }

        /* Integrate only active particles and save next positions/orientations in parts[i].next.
         * Non-active particles will usually correspond to one of the particles of a bound pair of particles */
        for (int i = 0; i < parts.size(); i++) {
            if (parts[i].isActive()) {
                /* Choose basic integration depending if particle MSM is
                 * active (this corresponds only to the MSM in the unbound state) */
                if (parts[i].activeMSM) {
                    integrateOneMS(i, parts, dt);
                } else {
                    integrateOne(i, parts, dt);
                }
            }
        }

        /* Advance global time and in event manager (to make events happen). Useful to draw a timeline to
         * understand order of events. */
        clock += dt;
        eventMgr.advanceTime(dt);

        // Remove unrealized previous events (see function for detailed description).
        removeUnrealizedEvents(parts);

        // Compute transitions to bound states (from unbound states) and add them to the event manager.
        computeTransitionsFromTransitionStates(parts);

        // Compute transitions from bound states (to unbound or other bound states) and add them to event manager.
        computeTransitionsFromBoundStates(parts);

        // Check for events in event manager that should happen during this time step [t,t+dt) and make them happen.
        applyEvents(parts);

        // Enforce boundary and set new positions into parts[i].nextPosition (only if particle is active).
        enforceBoundary(parts);

        /* Update positions and orientations (sets calculated next position/orientation
         * calculated by integrator and boundary as current position/orientation). Note states
         * are modified directly and don't need to be updated. */
        updatePositionOrientation(parts);

    }


}