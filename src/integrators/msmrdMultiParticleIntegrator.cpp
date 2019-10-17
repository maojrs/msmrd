//
// Created by maojrs on 10/7/19.
//

#include "integrators/msmrdMultiParticleIntegrator.hpp"

namespace msmrd {
    /**
     * Implementation of multiparticle MSM/RD integrator. Based mainly on msmrdintegratorDiscrete.
     * Uses same constructor as msmrdIntegratorDiscrete, with the followin parameters
     * @param dt time step
     * @param seed random generator seed (Note seed = -1 corresponds to random device)
     * @param rotation boolean to indicate if rotational degrees of freedom should be integrated
     */


    /* Computes possible transitions from transition states to bound states or other transition states from
     * particles sufficiently close to each other (in transition states) and saves them in the event manager.
     * Used by integrate function. */
    void msmrdMultiParticleIntegrator::computeTransitionsFromTransitionStates(std::vector<particleMS> &parts) {
        int currentTransitionState;
        double transitionTime;
        int nextState;
        int index0 = markovModel.getMaxNumberBoundStates();
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
    void msmrdMultiParticleIntegrator::computeTransitionsFromBoundStates(std::vector<particleMS> &parts) {
        double transitionTime;
        int nextState;
        std::tuple<double, int> transition;
        int index0 = markovModel.getMaxNumberBoundStates();
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
                        transition = markovModel.calculateTransition(parts[i].boundStates[boundParticleIndex]);
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


    // CURRENTLY WORKING HERE...

    /* Applies transition to bound state if binding pocket is not taken by another binding.
     * Makes particles with indexes iIndex and jIndex in the particle list transition to a bound state. Note
     * always iIndex < jIndex should hold. Also particle with smaller index is the one that remains active
     * to model position/orientation of bound complex. */
    void msmrdMultiParticleIntegrator::transition2BoundState(std::vector<particleMS> &parts, int iIndex,
                                                             int jIndex, int endState) {

        /* Establish pair connection in particle class and deactivate particle with larger index ( only
         * particle with smaller index remains active to represent the movement of the bound particle) */
        parts[iIndex].boundList.push_back(jIndex);
        parts[jIndex].boundList.push_back(iIndex);
        parts[jIndex].deactivate();
        // Set bound state for particle
        parts[iIndex].boundStates.push_back(endState);
        parts[jIndex].boundStates.push_back(endState);
        // Set diffusion coefficients in bound state (note states start counting from 1, not zero)
        int MSMindex = markovModel.getMSMindex(endState);
        parts[iIndex].setDs(markovModel.Dlist[MSMindex], markovModel.Drotlist[MSMindex]);
        // Average bound particle position and orientation (save on particle with smaller index).
        if (rotation) {
            parts[iIndex].nextOrientation = parts[iIndex].orientation;
            //parts[jIndex].nextOrientation = USE RELATIVE ORIENTATION HERE!!!!
        }
        parts[iIndex].nextPosition = parts[iIndex].position;
        // Assign constant distant position to inactive particle ( could be useful for visualization).
        parts[jIndex].nextPosition = parts[jIndex].position;
    }
}