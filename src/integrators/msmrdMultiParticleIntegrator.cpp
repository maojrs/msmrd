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

        // Add particle complex to particleComplexes vector.
        addComplex(parts, iIndex, jIndex, endState);

        // CURRENTLY WORKING HERE...

        // Set average positions and orientations of particles

        // Average bound particle position and orientation (save on particle with smaller index).
        if (rotation) {
            parts[iIndex].nextOrientation = parts[iIndex].orientation;
            //parts[jIndex].nextOrientation = USE RELATIVE ORIENTATION HERE!!!!
        }
        parts[iIndex].nextPosition = parts[iIndex].position;
        // Assign constant distant position to inactive particle ( could be useful for visualization).
        parts[jIndex].nextPosition = parts[jIndex].position;
    }





    /* Add particle complex into vector particleComplexes if particles bounded. If complex doesn't exist,
     * it creates it. If both particles belong to different complexes, it merges them. It also updates
     * the complexIndex of each particle, see particle.hpp.*/
    void msmrdMultiParticleIntegrator::addComplex(std::vector<particleMS> &parts, int iIndex,
                                                  int jIndex, int endState) {
        particleMS &iPart = parts[iIndex];
        particleMS &jPart = parts[jIndex];
        // If neither particle belongs to a complex, create one
        if (iPart.complexIndex == -1 and jPart.complexIndex == -1) {
            std::tuple<int,int> pairIndices = std::make_tuple(iIndex, jIndex);
            std::map<std::tuple<int,int>, int> boundPairsDictionary = {{pairIndices, endState}};
            particleComplex pComplex = particleComplex(iPart.position, iPart.orientation, boundPairsDictionary);
            particleComplexes.push_back(pComplex);
            // Set new particle complex indices.
            iPart.complexIndex = particleComplexes.size() - 1;
            jPart.complexIndex = particleComplexes.size() - 1;
        }
        // If one of the two doesn't belong to a complex, join the solo particle into the complex.
        else if (iPart.complexIndex * jPart.complexIndex < 0) {
            int complexIndex = std::min(iPart.complexIndex, jPart.complexIndex);
            std::tuple<int,int> pairIndices = std::make_tuple(iIndex, jIndex);
            particleComplexes[complexIndex].boundPairsDictionary.insert (
                    std::pair<std::tuple<int,int>, int>(pairIndices, endState) );
            // Set new particle complex indices.
            iPart.complexIndex = complexIndex;
            jPart.complexIndex = complexIndex;
        }
        // If both belong to a complex, join the complexes together (keep the one with lower index)
        else {
            // Add new bound pair to particle complex
            int iComplexIndex = std::min(iPart.complexIndex, jPart.complexIndex);
            int jComplexIndex = std::max(iPart.complexIndex, jPart.complexIndex);
            std::tuple<int,int> pairIndices = std::make_tuple(iIndex, jIndex);
            particleComplexes[iComplexIndex].boundPairsDictionary.insert (
                    std::pair<std::tuple<int,int>, int>(pairIndices, endState) );
            // Join complexes and flag complex with alrger index to be deleted.
            particleComplexes[iComplexIndex].joinParticleComplex(particleComplexes[jComplexIndex]);
            particleComplexes[jComplexIndex].active = false;
            // Set new particle complex indices to the one of the lower index.
            iPart.complexIndex = iComplexIndex;
            jPart.complexIndex = iComplexIndex;
        }
    };


    // @TODO WRITE TEST ROUTINE FOR THE FUNCTION BELOW: THIS FUNCTION WILL NEED TO BE ADDED IN THE INTEGRATOR FUNCTION

    /* Deletes inactive complexes in particle complex vector, and updates indexes in particle list. Doesn't
     * need to do at every time step, but every now and then to make free up memory. */
    void msmrdMultiParticleIntegrator::updateParticleComplexesVector(std::vector<particleMS> &parts) {
        for (size_t i = particleComplexes.size(); i--;) {
            if (particleComplexes[i].active == false) {
                // Erase particle complex
                particleComplexes.erase(particleComplexes.begin() + i);
                // Readjust the indexes of particles
                for (int j = 0; j < parts.size(); j++) {
                    // this should never happen
                    if (parts[j].complexIndex == i) {
                        parts[j].complexIndex = -1;
                    }
                    // move all the indices larger than i by -1 since we will delete one element
                    if (parts[j].complexIndex > i) {
                        parts[j].complexIndex--;
                    }
                }
            }
        }
    };



}