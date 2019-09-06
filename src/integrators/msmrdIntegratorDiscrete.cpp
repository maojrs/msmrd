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
                //relativeOrientation = part1.nextOrientation.conj() * part2.nextOrientation;
                relativeOrientation = part2.nextOrientation * part1.nextOrientation.conj();
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
                    if (previousEvent.eventType == "empty") {
                        // returns -1 if |relativePosition| > radialBounds[1]
                        currentTransitionState = computeCurrentTransitionState(parts[i], parts[j]);
                    } else if (previousEvent.eventType == "inTransition") {
                        //previous endState is current starting state
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
        int MSMindex = markovModel.getMSMindex(endState);
        parts[iIndex].setDs(markovModel.Dlist[MSMindex], markovModel.Drotlist[MSMindex]);
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

        /* Calculate new relative positions and orientations by sampling either:
         * nonuniformly on sphere section weighted on r close to the center,
         * uniformly in spherical section, uniformly in outer shell section,
         * uniformly in inner shell section or custom */
        //auto rr = randg.uniformRange(radialBounds[0],radialBounds[1]);
        //auto relPosition = rr * randg.uniformSphereSection(phiInterval, thetaInterval);
        //auto relPosition = randg.uniformShellSection(radialBounds, phiInterval, thetaInterval);
        //auto relPosition = radialBounds[1] * randg.uniformSphereSection(phiInterval, thetaInterval);
        auto relPosition = radialBounds[0] * randg.uniformSphereSection(phiInterval, thetaInterval);
        //double rr  = radialBounds[0] + 0.2*(radialBounds[1] - radialBounds[0]);
        //auto relPosition = rr * randg.uniformSphereSection(phiInterval, thetaInterval);


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
        /* Change eventType label to indicate to computeTransitionsFromTransitionStates function that
         * a new event needs to be calculated, using previousEvent.endState as the initial state */
        eventMgr.setEventType("inTransition", iIndex, jIndex);
    }


    /* Removes unrealized events where unbound particles drifted a distance apart beyond the upper radial bound,
     * or when zero rates yielded infinite values. */
    template<>
    void msmrdIntegratorDiscrete<ctmsm>::removeUnrealizedEvents(std::vector<particleMS> &parts) {
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
    void msmrdIntegratorDiscrete<ctmsm>::applyEvents(std::vector<particleMS> &parts) {
        std::list<std::map<std::string, decltype(eventMgr.emptyEvent)>::const_iterator> iteratorList;
        // Loop over dictionary using iterators
        for (auto it = eventMgr.eventDictionary.cbegin(); it != eventMgr.eventDictionary.cend(); it++) {
            auto transitionTime = it->second.waitTime;
            /* Only apply events that should happen in during this
             * timestep (they correspond to transitionTime < 0) */
            if (transitionTime <= 0) {
                // Load event data
                auto iIndex = it->second.part1Index;
                auto jIndex = it->second.part2Index;
                auto endState = it->second.endState;
                auto eventType = it->second.eventType;
                // Make event happen (depending on event type) and remove event once it has happened
                if (eventType == "binding") {
                    transition2BoundState(parts, iIndex, jIndex, endState);
                    iteratorList.push_back(it);
                } else if (eventType == "unbinding") {
                    transition2UnboundState(parts, iIndex, jIndex, endState);
                    iteratorList.push_back(it);
                } else if (eventType == "bound2bound") {
                    transitionBetweenBoundStates(parts, iIndex, jIndex, endState);
                    iteratorList.push_back(it);
                } else if (eventType == "transition2transition") {
                    /* Note event is not removed until a new event is computed later in
                     * the computeTransitionsFromTransitionStates routine */
                    transitionBetweenTransitionStates(iIndex, jIndex);
                }
            }
        }
        // Erase events that just were applied
        for (auto it : iteratorList) {
            eventMgr.eventDictionary.erase(it);
        }
    }

    /* Integrates translational and rotation diffusion for on dt step. (Calls the integrateOne fucntions
     * of the parent classes) */
    template<>
    void msmrdIntegratorDiscrete<ctmsm>::integrateDiffusion(std::vector<particleMS> &parts, double dt) {
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

        /* NOTE: the ordering of the following routines is veryy important, draw a timeline if necessary.*/

        // Compute future transitions to bound states (from unbound states) and add them to the event manager.
        computeTransitionsFromTransitionStates(parts);

        // Compute future transitions from bound states (to unbound or other bound states); add them to event manager.
        computeTransitionsFromBoundStates(parts);

        // Integrates diffusion for one time step.
        integrateDiffusion(parts, dt);

        // Remove unrealized previous events (see function for detailed description).
        removeUnrealizedEvents(parts);

        /* Advance global time and in event manager (to make events happen). Useful to draw a timeline to
         * understand order of events. */
        clock += dt;
        eventMgr.advanceTime(dt);

        /* Check for events in event manager that should happen during this time step [t,t+dt) and
         * make them happen. Note if the lag-time is a multiple of dt (n*dt), which is likely the case,
         * events will happen at end of the timestep exactly */
        applyEvents(parts);

        // Enforce boundary and set new positions into parts[i].nextPosition (only if particle is active).
        enforceBoundary(parts);

        /* Update positions and orientations (sets calculated next position/orientation
         * calculated by integrator and boundary as current position/orientation). Note states
         * are modified directly and don't need to be updated. */
        updatePositionOrientation(parts);


    }


}