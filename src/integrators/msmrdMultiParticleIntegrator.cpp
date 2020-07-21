//
// Created by maojrs on 10/7/19.
//

#include "integrators/msmrdMultiParticleIntegrator.hpp"

namespace msmrd {
    /**
     * Implementation of multiparticle MSM/RD integrator. Based mainly on msmrdintegrator.
     * Uses same constructor as msmrdIntegrator, with the following parameters
     * @param dt time step
     * @param seed random generator seed (Note seed = -1 corresponds to random device)
     * @param rotation boolean to indicate if rotational degrees of freedom should be integrated
     */

    /* Computes possible transitions from transition states to bound states or other transition states from
     * particles sufficiently close to each other (in transition states) and saves them in the event manager.
     * Used by integrate function. */
    template<>
    void msmrdMultiParticleIntegrator<ctmsm>::computeTransitionsFromTransitionStates(std::vector<particle> &parts) {
        int currentTransitionState;
        double transitionTime;
        int nextState;
        int index0 = msmrdMSM.getMaxNumberBoundStates();
        // Loop over all pairs of particles without repetition with i < j
        for (int i = 0; i < parts.size(); i++) {
            for (int j = i + 1; j < parts.size(); j++) {
                /* Only compute transitions if both particles have at least one bound site free (bound to one or zero
                 * other particles) */
                if (parts[i].boundList.size() < 2 and parts[j].boundList.size() < 2) {
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
                        auto transition = msmrdMSM.calculateTransition(currentTransitionState);
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


    /* NOT USED IN FIRST IMPLEMENTTION .
     * Computes possible transitions from bound states to other bound states or unbound states (transition states)
     * and saves them in the event manager. Used by integrate function. */
    template<>
    void msmrdMultiParticleIntegrator<ctmsm>::computeTransitionsFromBoundStates(std::vector<particle> &parts) {
        double transitionTime;
        int nextState;
        std::tuple<double, int> transition;
        int index0 = msmrdMSM.getMaxNumberBoundStates();
        for (int i = 0; i < parts.size(); i++) {
            // Only compute transition if particle is bound to another particle.
            // If particle[i] is bound, then parts[i].boundList is not empty.
            for (auto boundParticleIndex : parts[i].boundList) {
                // Bound pairs are only to be counted once, then:
                if (parts[i].boundList[boundParticleIndex] > i) {
                    /* Only compute transition if particles switched into a given bound state for
                     * the first time, i.e. empty event */
                    auto previousEvent = eventMgr.getEvent(i, parts[i].boundList[boundParticleIndex]);
                    if (previousEvent.eventType == "empty") {
                        transition = msmrdMSM.calculateTransition(parts[i].boundStates[boundParticleIndex]);
                        transitionTime = std::get<0>(transition);
                        nextState = std::get<1>(transition);
                        // Distinguish between events bound to bound transition and unbinding events
                        if (nextState <= index0) {
                            eventMgr.addEvent(transitionTime, i, parts[i].boundList[boundParticleIndex],
                                              parts[i].boundStates[boundParticleIndex], nextState, "bound2bound");
                        } else {
                            eventMgr.addEvent(transitionTime, i, parts[i].boundList[boundParticleIndex],
                                              parts[i].boundStates[boundParticleIndex], nextState, "unbinding");
                        }
                    }
                }
            }
        }
    }


    /* Applies transition to bound state if binding pocket is not taken by another binding.
     * Makes particles with indexes iIndex and jIndex in the particle list transition to a bound state. Note
     * always iIndex < jIndex should hold. Also both particles are desactivated since their position will be
     * tracked with a particle compound. */
    template<>
    void msmrdMultiParticleIntegrator<ctmsm>::transition2BoundState(std::vector<particle> &parts, int iIndex,
                                                             int jIndex, int endState) {
        /* Establish pair connection in particle class and deactivate both particles */
        parts[iIndex].boundList.push_back(jIndex);
        parts[jIndex].boundList.push_back(iIndex);
        parts[iIndex].deactivate();
        parts[jIndex].deactivate();
        // Set bound state for particle
        parts[iIndex].boundStates.push_back(endState);
        parts[jIndex].boundStates.push_back(endState);

        // Add particle complex to particleComplexes vector.
        int mainCompoundSize = addCompound(parts, iIndex, jIndex, endState);

        // Set average position and orientation of particle compound
        setCompoundPositionOrientation(parts, iIndex, jIndex, mainCompoundSize);

        /* Set diffusion coefficients of bound particle compound (note states start counting from 1, not zero).
         * In general this should be a more complicated when binding larger complexes. */
        int MSMindex = msmrdMSM.getMSMindex(endState);
        if (particleCompounds[parts[iIndex].compoundIndex].compoundSize == 2) {
            particleCompounds[parts[iIndex].compoundIndex].setDs(msmrdMSM.Dlist[MSMindex],
                                                                 msmrdMSM.Drotlist[MSMindex]);
        }
        /* TODO: implement diffusion coeffcients assignments for larger complexes (> 2 particles).
         * Not too relevant in first implementation since all diffusion coefficients are constats,
         * so left for future implemenation. */
    }


    /* NOT ALLOWING FOR UNBINDING EVENTS IN FIRST IMPLEMENTATION.
     * CURRENTLY WORKING ON FUNCTION BELOW, ENCOUNTERED ISSUE IN HOW PARTICLECOMPOUND SAVES TOPOLOGIES, NEED TO
     * ADDRESS THAT FIRST, SEE particle.hpp */
    template<>
    void msmrdMultiParticleIntegrator<ctmsm>::transition2UnboundState(std::vector<particle> &parts, int iIndex,
                                                               int jIndex, int endStateAlt) {

        // Check if particles unbinding belong to particleCompounds
        parts[iIndex].compoundIndex;
        parts[iIndex].compoundIndex;

        // Redefine endstate indexing, so it is understood by the partition/discretization.
        int index0 = msmrdMSM.getMaxNumberBoundStates();
        int endState = endStateAlt - index0;

        // Calculates and sets next unbound states (of the unbound MSM). If no MSM, defaults to zero.
        setRandomUnboundState(parts, iIndex);
        setRandomUnboundState(parts, jIndex);

        // Extract relative position and orientation from partition and endstate
        auto relativePositionOrientation = getRelativePositionOrientation(endState);
        auto relPosition = std::get<0>(relativePositionOrientation);

        // Set next orientations based on the relative ones (parts[iIndex] keeps track of bound particle orientation)
        if (rotation) {
            auto relOrientation = std::get<1>(relativePositionOrientation);
            parts[iIndex].nextOrientation = 1.0 * parts[iIndex].orientation;
            parts[jIndex].nextOrientation = relOrientation * parts[iIndex].nextOrientation;
        }

        // Set next positions based on the relative ones (parts[iIndex] keeps track of bound particle position)
        parts[iIndex].nextPosition = parts[iIndex].position - 0.5*relPosition;
        parts[jIndex].nextPosition = parts[iIndex].nextPosition + relPosition;
    }


    /* Main integrate function */
    template<>
    void msmrdMultiParticleIntegrator<ctmsm>::integrate(std::vector<particle> &parts) {

        /* Calculate forces and torques and save them into forceField and torqueField. For the MSM/RD this will
         * in general be zero, so only needs to be run once. However, in some case they might be activated */
        if (firstrun or pairPotentialActive or externalPotentialActive) {
            calculateForceTorqueFields<particle>(parts);
            firstrun = false;
        }

        /* NOTE: the ordering of the following routines is very important, draw a timeline if necessary.*/

        // Compute future transitions to bound states (from unbound states) and add them to the event manager.
        computeTransitionsFromTransitionStates(parts);

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


        // Output eventlog (useful for debugging)
        if (recordEventLog) {
            auto timeIteration = static_cast<int>(clock / dt);
            eventMgr.write2EventLog(timeIteration);
        }

    }



    /**
     * Additional functions exclusive to multi-particle MSM/RD below
     */

    // @TODO WRITE TEST ROUTINE FOR THE FUNCTION BELOW: THIS FUNCTION WILL NEED TO BE ADDED IN THE INTEGRATOR FUNCTION
    template<>
    void msmrdMultiParticleIntegrator<ctmsm>::integrateDiffusionCompounds(double dt0,
                                                                                std::vector<particle> &parts){
        for (auto particleCompound : particleCompounds) {
            // Calculte change in poisition
            vec3<double> dr;
            dr =  std::sqrt(2 * dt0 * particleCompound.D) * randg.normal3D(0, 1);
            particleCompound.position += dr;
            // Calculate change in orientation
            vec3<double> dphi;
            quaternion<double> dquat;
            dphi = std::sqrt(2 * dt0 * particleCompound.Drot) * randg.normal3D(0, 1);
            dquat = msmrdtools::axisangle2quaternion(dphi);
            particleCompound.orientation = dquat * particleCompound.orientation;
            // Make a copy of the boundsPair dictionary
            std::map<std::tuple<int,int>, int> boundPairsDictionaryCopy(particleCompound.boundPairsDictionary);
            // Update position and orientation of reference particle
            setParticlePositionOrientation(parts, particleCompound, dr,  dquat, particleCompound.referenceIndex);
            /* Update position and orientation of all other particles in compound (depends on
             * specific implementation, here based on the dimer example with one angular binding). Loops many times
             * over copy of boundsPairsDictionary, starting by the reference particle and building all the bindings
             * until there is no binding left to incorporate. */
            int tempIndex1 = 1 * particleCompound.referenceIndex;
            int tempIndex2 = 1 * particleCompound.referenceIndex;
            while (boundPairsDictionaryCopy.size() > 0) {
                for (auto it = boundPairsDictionaryCopy.cbegin(); it != boundPairsDictionaryCopy.cend(); it++) {
                    auto index1 = std::get<0>(it->first);
                    auto index2 = std::get<1>(it->first);
                    auto state = it->second;
                    if (index1 == tempIndex1) {
                        // Assign positions and orientations to particles involved
                        tempIndex1 = index2;
                    }
                    if (index2 == tempIndex2) {
                        // Assign positions and orientations to particles involved
                        tempIndex2 = index1;
                    }
                    // Remove element from dictionary and break loop
                    boundPairsDictionaryCopy.erase(it);
                    break;
                }
            }
        }
    }







}