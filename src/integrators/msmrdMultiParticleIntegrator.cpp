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
    void msmrdMultiParticleIntegrator::computeTransitionsFromTransitionStates(std::vector<particle> &parts) {
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


    /* Computes possible transitions from bound states to other bound states or unbound states (transition states)
     * and saves them in the event manager. Used by integrate function. */
    void msmrdMultiParticleIntegrator::computeTransitionsFromBoundStates(std::vector<particle> &parts) {
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
    void msmrdMultiParticleIntegrator::transition2BoundState(std::vector<particle> &parts, int iIndex,
                                                             int jIndex, int endState) {
        /* Establish pair connection in particle class and deactivate particle with larger index ( only
         * particle with smaller index remains active to represent the movement of the bound particle) */
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


    /* CURRENTLY WORKING HERE BELOW, ENCOUNTERED ISSUE IN HOW PARTICLECOMPOUND SAVES TOPOLOGIES, NEED TO
     * ADDRESS THAT FIRST, SEE particle.hpp */

    void msmrdMultiParticleIntegrator::transition2UnboundState(std::vector<particle> &parts, int iIndex,
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



    /**
     * Functions exclusive to multi-particle MSM/RD below
     */


    /* Add particle complex into vector particleComplexes if particles bounded. If complex doesn't exist,
     * it creates it. If both particles belong to different complexes, it merges them. It also updates
     * the compoundIndex of each particle, see particle.hpp. It also returns an integer corresponding to
     * the size of the main compound in the case of two compounds joining, or simply one otherwise. */
    int msmrdMultiParticleIntegrator::addCompound(std::vector<particle> &parts, int iIndex,
                                                  int jIndex, int endState) {
        int mainCompoundSize = 1;
        particle &iPart = parts[iIndex];
        particle &jPart = parts[jIndex];
        // If neither particle belongs to a complex, create one
        if (iPart.compoundIndex == -1 and jPart.compoundIndex == -1) {
            std::tuple<int,int> pairIndices = std::make_tuple(iIndex, jIndex);
            std::map<std::tuple<int,int>, int> boundPairsDictionary = {{pairIndices, endState}};
            particleCompound pComplex = particleCompound(boundPairsDictionary);
            particleCompounds.push_back(pComplex);
            // Set new particle complex indices.
            iPart.compoundIndex = static_cast<int>(particleCompounds.size() - 1);
            jPart.compoundIndex = static_cast<int>(particleCompounds.size() - 1);
        }
        // If one of the two doesn't belong to a complex, join the solo particle into the complex.
        else if (iPart.compoundIndex * jPart.compoundIndex < 0) {
            // Generate new binding description
            int compoundIndex = std::min(iPart.compoundIndex, jPart.compoundIndex);
            std::tuple<int,int> pairIndices = std::make_tuple(iIndex, jIndex);
            //Insert new binding decription into complex
            particleCompounds[compoundIndex].boundPairsDictionary.insert (
                    std::pair<std::tuple<int,int>, int>(pairIndices, endState) );
            // Extract main compound size and increase complex size by one.
            mainCompoundSize = particleCompounds[compoundIndex].compoundSize;
            particleCompounds[compoundIndex].compoundSize ++;
            // Set new particle complex indices.
            iPart.compoundIndex = compoundIndex;
            jPart.compoundIndex = compoundIndex;
        }
        // If both belong to a complex, join the complexes together (keep the one with lower index)
        else {
            // Add new bound pair to particle complex
            int iCompoundIndex = std::min(iPart.compoundIndex, jPart.compoundIndex);
            int jCompoundIndex = std::max(iPart.compoundIndex, jPart.compoundIndex);
            std::tuple<int,int> pairIndices = std::make_tuple(iIndex, jIndex);
            particleCompounds[iCompoundIndex].boundPairsDictionary.insert (
                    std::pair<std::tuple<int,int>, int>(pairIndices, endState) );
            // Join complexes and flag complex with alrger index to be deleted.
            particleCompounds[iCompoundIndex].joinParticleCompound(particleCompounds[jCompoundIndex]);
            particleCompounds[jCompoundIndex].active = false;
            // Extract main compound size, increase complex size by one and make zero the empty one
            mainCompoundSize = particleCompounds[iCompoundIndex].compoundSize;
            particleCompounds[iCompoundIndex].compoundSize ++;
            particleCompounds[jCompoundIndex].compoundSize = 0;
            // Set new particle complex indices to the one of the lower index.
            iPart.compoundIndex = iCompoundIndex;
            jPart.compoundIndex = iCompoundIndex;
        }
        return mainCompoundSize;
    };


    /* Sets compound position to the average, its orientation is always set up initially at one and only changed
     * by diffusion. It also sets the position reference and orientation reference for a newly created compound. */
    void msmrdMultiParticleIntegrator::setCompoundPositionOrientation(std::vector<particle> &parts,
                                                                      int iIndex, int jIndex, int mainCompoundSize) {
        int compoundIndex = parts[iIndex].compoundIndex;
        int compoundSize = particleCompounds[compoundIndex].compoundSize;
        if (compoundSize == 2) {
            // Set average position
            particleCompounds[compoundIndex].position = 0.5*(parts[iIndex].position + parts[jIndex].position);
            vec3<double> diffVec = parts[jIndex].position - parts[iIndex].position;
            particleCompounds[compoundIndex].positionReference = parts[iIndex].position - 0.5*diffVec/diffVec.norm();
            particleCompounds[compoundIndex].orientationReference = 1.0*parts[iIndex].orientation;
        } else {
            particleCompounds[compoundIndex].position = (mainCompoundSize * particleCompounds[compoundIndex].position +
                                                         (compoundSize - mainCompoundSize) * parts[jIndex].position)/compoundSize;
        }
    };


//    /* Checks ...*/
//    std::tuple<bool, bool> msmrdMultiParticleIntegrator::checkUnbindingCompounds(std::vector<particle> &parts,
//                                                                               int iIndex, int jIndex) {
//        int cIndex = parts[iIndex].compoundIndex;
//
//        for (const auto &entry : particleCompounds[cIndex].boundPairsDictionary) {
//            auto key = entry.first;
//            auto value = entry.second
//        }
//
//
//    };



    // @TODO WRITE TEST ROUTINE FOR THE FUNCTION BELOW: THIS FUNCTION WILL NEED TO BE ADDED IN THE INTEGRATOR FUNCTION

    /* Deletes inactive complexes in particle complex vector, and updates indexes in particle list. Doesn't
     * need to do at every time step, but every now and then to free up memory. */
    void msmrdMultiParticleIntegrator::updateParticleComplexesVector(std::vector<particle> &parts) {
        for (size_t i = particleCompounds.size(); i--;) {
            if (particleCompounds[i].active == false) {
                // Erase particle complex
                particleCompounds.erase(particleCompounds.begin() + i);
                // Readjust the indexes of particles
                for (int j = 0; j < parts.size(); j++) {
                    // this should never happen
                    if (parts[j].compoundIndex == i) {
                        parts[j].compoundIndex = -1;
                    }
                    // move all the indices larger than i by -1 since we will delete one element
                    if (parts[j].compoundIndex > i) {
                        parts[j].compoundIndex--;
                    }
                }
            }
        }
    };



}