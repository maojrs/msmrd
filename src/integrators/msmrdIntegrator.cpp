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
        // Loop over all pairs of particles without repetition with i < j
        for (int i = 0; i < parts.size(); i++) {
            for (int j = i + 1; j < parts.size(); j++) {
                // Only compute transitions if particles are in unbound state.
                if (parts[i].boundTo == -1 and parts[j].boundTo == -1) {
                    /* Only compute new transition if particles drifted into transition region for
                     * the first time, i.e. empty event and relativeDistance < radialBounds[1] */
                    auto previousEvent = eventMgr.getEvent(i, j);
                    if (previousEvent.eventType == "empty") {
                        // Need special function to calculate relative position, in case we have a periodic boundary.
                        relativePosition = calculateRelativePosition(parts[i].nextPosition, parts[j].nextPosition);
                        if (relativePosition.norm() < radialBounds[1]) {
                            if (rotation) {
                                //relativeOrientation = parts[i].nextOrientation.conj() * parts[j].nextOrientation;
                                relativeOrientation = parts[j].nextOrientation * parts[i].nextOrientation.conj();
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
                            eventMgr.addEvent(transitionTime, i, j, currentTransitionState, nextState, "binding");
                        }
                    }
                }
            }
        }
    }

    /* Computes possible transitions from bound states to other bound states or unbound states (transition states)
     * and saves them in the event manager. Used by integrate function. */
    template<>
    void msmrdIntegrator<ctmsm>::computeTransitionsFromBoundStates(std::vector<particleMS> &parts) {
        double transitionTime;
        int nextState;
        std::tuple<double, int> transition;
        for (int i = 0; i < parts.size(); i++) {
            // Only compute transition if particle is bound to another particle.
            // If particle[i] is bound, boundTo > 0; if bound pairs are only to be counted once, then boundTo > i
            if (parts[i].boundTo > i) {
                /* Only compute transition if particles switched into a given bound state for
                 * the first time, i.e. empty event */
                auto previousEvent = eventMgr.getEvent(i, parts[i].boundTo);
                if (previousEvent.eventType == "empty") {
                    transition = markovModel.computeTransitionFromBoundState(parts[i].state);
                    transitionTime = std::get<0>(transition);
                    nextState = std::get<1>(transition);
                    // Distinguish between events bound to bound transition) and unbinding events
                    if (nextState <= markovModel.getMaxNumberBoundStates()) {
                        eventMgr.addEvent(transitionTime, i, parts[i].boundTo,
                                          parts[i].state, nextState, "bound2boundTransition");
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
        // Redefine endstate indexing, so it is understood by the partition/discretization.
        int maxNumberBoundStates = markovModel.getMaxNumberBoundStates();
        endState = endState - maxNumberBoundStates;
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

        /* Calculate new relative positions and orientations by sampling either:
         * nonuniformly on sphere section weighted on r close to the center,
         * uniformly in spherical section, uniformly in outer shell section or
         * uniformly in inner shell section. */
        //auto rr = randg.uniformRange(radialBounds[0],radialBounds[1]);
        //auto relPosition = rr * randg.uniformSphereSection(phiInterval, thetaInterval);
        //auto relPosition = randg.uniformShellSection(radialBounds, phiInterval, thetaInterval);
        //auto relPosition = radialBounds[1] * randg.uniformSphereSection(phiInterval, thetaInterval);
        auto relPosition = radialBounds[0] * randg.uniformSphereSection(phiInterval, thetaInterval);

        // Set next positions and orientations based on the relative ones (parts[iIndex] keeps track of position)
        parts[iIndex].nextPosition = parts[iIndex].position - 0.5*relPosition;
        parts[jIndex].nextPosition = parts[iIndex].nextPosition + relPosition;
    }


    /* Transitions of already bound particles with indexes iIndex and jIndex in the particle list
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


    /* Removes unrealized events where unbound particles drifted a distance apart beyond the upper radial bound,
     * or when zero rates yielded infinite values. */
    template<>
    void msmrdIntegrator<ctmsm>::removeUnrealizedEvents(std::vector<particleMS> &parts) {
        std::list<std::map<std::string, decltype(eventMgr.emptyEvent)>::const_iterator> iteratorList;
        vec3<double> relativePosition;
        // Loop over all events to flag events that should be erased (loop over iterator)
        for (auto it = eventMgr.eventDictionary.cbegin(); it != eventMgr.eventDictionary.cend(); it++) {
            auto transitionTime = it->second.waitTime;
            auto iIndex = it->second.part1Index;
            auto jIndex = it->second.part2Index;
            // Flag event to be removed if transition time is infinity
            if (std::isinf(transitionTime)) {
                iteratorList.push_back(it);
            }
            // If particles in unbound state and relative position larger than cutOff, flag event to be removed.
            if (parts[iIndex].boundTo == -1 and parts[jIndex].boundTo == -1) {
                relativePosition = calculateRelativePosition(parts[iIndex].nextPosition, parts[jIndex].nextPosition);
                // Remove event if particles drifted apart
                if (relativePosition.norm() >= radialBounds[1]) {
                    iteratorList.push_back(it);
                }
            }
        }
        // Erase events flagged to be erased
        for (auto it : iteratorList) {
            eventMgr.eventDictionary.erase(it);
        }
    }


    /* Apply events in event manager that should happen during the current time step. */
    template<>
    void msmrdIntegrator<ctmsm>::applyEvents(std::vector<particleMS> &parts) {
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
                // Make event happen (depending on event type)
                if (eventType == "binding") {
                    transition2BoundState(parts, iIndex, jIndex, endState);
                } else if (eventType == "unbinding") {
                    transition2UnboundState(parts, iIndex, jIndex, endState);
                } else if (eventType == "bound2boundTransition") {
                    transitionBetweenBoundStates(parts, iIndex, jIndex, endState);
                }
                // Remove event from event list once it has happened
                eventMgr.removeEvent(iIndex, jIndex);
            }
        }
    }


    /* Main integrate function */
    template<>
    void msmrdIntegrator<ctmsm>::integrate(std::vector<particleMS> &parts) {

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
        computeTransitions2BoundStates(parts);

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