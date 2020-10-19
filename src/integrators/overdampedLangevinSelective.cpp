//
// Created by maojrs on 10/15/20.
//


#include "integrators/overdampedLangevinSelective.hpp"

namespace msmrd {

    overdampedLangevinSelective::overdampedLangevinSelective(double dt, long seed, std::string particlesbodytype) :
            overdampedLangevin(dt, seed, particlesbodytype) {
        patchyTraj = std::make_shared<discreteTrajectory<4>> (patchyDimerTrajectory2(2,1));
        patchyTraj->setTolerances(0.50, 0.5*2*M_PI);
    }

    /* Integrate list of particles with selective patch selection, see setActivePatches fucntion */
    void overdampedLangevinSelective::integrate(std::vector<particle> &parts) {
        vec3<double> force;
        vec3<double> torque;

        // Calculate forces and torques and save them into forceField and torqueField
        calculateForceTorqueFields(parts);

        for (int i = 0; i < parts.size(); i++) {
            // Sets active patches to calculate force-field avoind triple bindings
            setActivePatches(parts);
            // Integrate and save next positions/orientations in parts[i].next***
            integrateOne(i, parts, dt);
        }

        // Enforce boundary; sets new positions into parts[i].nextPosition (only if particle is active)
        enforceBoundary(parts);

        // Update position/orientation based on parts[i].nextPosition, parts[i].nextOrientation
        updatePositionOrientation(parts);

        // Updates time
        clock += dt;
    }


    /* Check and set which patches should be active or inactive. Any patch involved in a binding
     * will be deactivated and only allowed to interact with the bound particle. */
    void overdampedLangevinSelective::setActivePatches(std::vector<particle> &parts) {
        // Set active patches to true and clear bound lists
        for (int i = 0; i < parts.size(); i++) {
            for (int j=0; j < parts[i].activePatchList.size(); j++) {
                parts[i].activePatchList[j] = -1;
            }
        }
         /* Deactivate patches and set bound lists depending if particles are bound by looping
          * over all possible pairs */
        for (int i = 0; i < parts.size() - 1; i++) {
            for (int j = i+1; j < parts.size(); j++) {
                // Only allow a new bound state if there is no previous bound state between
                auto state = patchyTraj->sampleDiscreteState(parts[i], parts[j]);
                if (state == 1) {
                    if (parts[i].activePatchList[0] < 0 and parts[j].activePatchList[1] < 0) {
                        parts[i].activePatchList[0] = j;
                        parts[j].activePatchList[1] = i;
                    }
                } else if (state == 2) {
                    if (parts[i].activePatchList[0] < 0 and parts[j].activePatchList[0] < 0) {
                        parts[i].activePatchList[0] = j;
                        parts[j].activePatchList[0] = i;
                    }
                } else if (state == 3) {
                    if (parts[i].activePatchList[1] < 0 and parts[j].activePatchList[0] < 0) {
                        parts[i].activePatchList[1] = j;
                        parts[j].activePatchList[0] = i;
                    }
                } else if (state == 4) {
                    if (parts[i].activePatchList[1] < 0 and parts[j].activePatchList[1] < 0) {
                        parts[i].activePatchList[1] = j;
                        parts[j].activePatchList[1] = i;
                    }
                }
            }
        }
    }

    /* Updates the vector of particle compounds. Whenever it is calles, it erases the vector
     * of particleCompounds and repopulates it depending on the current configuration. Note it
     * doesn't keep track of actual realtive positions nor orientations. It is only to track bindings. */
    void overdampedLangevinSelective::updateParticleCompounds(std::vector<particle> &parts) {
        auto dummyRelPosition = vec3<double>(0,0,0);
        particleCompounds.clear();
        for (int i = 0; i < parts.size(); i++) {
            parts[i].compoundIndex = -1;
        }
        for (int i = 0; i < parts.size() - 1; i++) {
            for (int j = i+1; j < parts.size(); j++) {
                // Only allow a new bound state if there is no previous bound state between
                int state = patchyTraj->sampleDiscreteState(parts[i], parts[j]);
                if (state == 1 or state == 2 or state == 3 or state == 4) {
                    // Create compound
                    if (parts[i].compoundIndex == -1 and parts[j].compoundIndex == -1) {
                        std::tuple<int, int> pairIndices = std::make_tuple(i, j);
                        std::map<std::tuple<int, int>, int> boundPairsDictionary = {{pairIndices, state}};
                        particleCompound pComplex = particleCompound(boundPairsDictionary);
                        pComplex.relativePositions.insert(std::pair<int, vec3<double>>(i, dummyRelPosition));
                        pComplex.relativePositions.insert(std::pair<int, vec3<double>>(j, dummyRelPosition));
                        particleCompounds.push_back(pComplex);
                        parts[i].compoundIndex = static_cast<int>(particleCompounds.size() - 1);
                        parts[j].compoundIndex = static_cast<int>(particleCompounds.size() - 1);
                    }
                    // Add particle to compound
                    else if (parts[i].compoundIndex == -1 or parts[j].compoundIndex == -1) {
                        if (parts[i].compoundIndex != -1){
                            auto compIndex = parts[i].compoundIndex;
                            std::tuple<int,int> pairIndices = std::make_tuple(i,j);
                            particleCompounds[compIndex].boundPairsDictionary.insert (
                                    std::pair<std::tuple<int,int>, int>(pairIndices, state));
                            particleCompounds[compIndex].relativePositions.insert(
                                    std::pair<int, vec3<double>>(j,1 * dummyRelPosition));
                            parts[j].compoundIndex = 1 * parts[i].compoundIndex;
                        }
                        else{
                            auto compIndex = parts[j].compoundIndex;
                            std::tuple<int,int> pairIndices = std::make_tuple(j,i);
                            particleCompounds[compIndex].boundPairsDictionary.insert (
                                    std::pair<std::tuple<int,int>, int>(pairIndices, state));
                            particleCompounds[compIndex].relativePositions.insert(
                                    std::pair<int, vec3<double>>(i,1 * dummyRelPosition));
                            parts[i].compoundIndex = 1 * parts[j].compoundIndex;
                        }
                    }
                    // Join compounds
                    else if (parts[i].compoundIndex != parts[j].compoundIndex) {
                        auto compIndex = parts[i].compoundIndex;
                        auto secondCompIndex = parts[j].compoundIndex;
                        std::tuple<int,int> pairIndices = std::make_tuple(i, j);
                        particleCompounds[compIndex].boundPairsDictionary.insert (
                                std::pair<std::tuple<int,int>, int>(pairIndices, state));
                        particleCompounds[compIndex].joinCompound(particleCompounds[secondCompIndex]);
                        particleCompounds[compIndex].relativePositions.insert(
                                particleCompounds[secondCompIndex].relativePositions.begin(),
                                particleCompounds[secondCompIndex].relativePositions.end());
                        particleCompounds[secondCompIndex].deactivateCompound();
                        for (int k=0; k < parts.size(); k++) {
                            if(parts[k].compoundIndex == parts[j].compoundIndex){
                                parts[k].compoundIndex = compIndex;
                            }
                        }
                    }
                    // Close compound
                    else{
                        std::tuple<int,int> pairIndices = std::make_tuple(i, j);
                        particleCompounds[parts[i].compoundIndex].boundPairsDictionary.insert (
                                std::pair<std::tuple<int,int>, int>(pairIndices, state) );
                    }
                }
            }
        }
    }

    /* Checks if there is any closed binding loop in any of the particle compounds. If so, it returns the size of
     * the loops found in a vector of integers, which size is the number of loops */
    std::vector<int> overdampedLangevinSelective::findClosedBindingLoops(std::vector<particle> &parts){
        // Updates particle compounds array to keep track of bindings.
        updateParticleCompounds(parts);
        // Check for ring loops in updated particle compounds vector.
        std::vector<int> boundLoops;
        for (auto &particleCompound : particleCompounds) {
            if (particleCompound.active) {
                auto compoundSize = particleCompound.getSizeOfCompound();
                auto numBindings = particleCompound.getNumberOfbindings();
                if (compoundSize == numBindings) {
                    boundLoops.push_back(numBindings);
                }
            }
        }
        return boundLoops;
    };


    int overdampedLangevinSelective::getCompoundSize(int compoundIndex) {
        if (compoundIndex < particleCompounds.size()) {
            return particleCompounds[compoundIndex].getSizeOfCompound();
        }
        else {
            return -1;
        }
    };



//    /* Check if ring molecule was formed. It returns 0 if it was not formed or it returns an
//     * indteger indicating the number of molecules involved in the ring.  */
//    int overdampedLangevinSelective::hasRingFormed(std::vector<particle> &parts) {
//        auto pentamerFormed = false;
//        std::vector<int> bindingsListPatch1(5, 0);
//        std::vector<int> bindingsListPatch2(5, 0);
//        // Count bindings in patches
//        for (int i = 0; i < parts.size() - 1; i++) {
//            for (int j = i + 1; j < parts.size(); j++) {
//                // Only allow a new bound state if there is no previous bound state between
//                auto state = patchyTraj->sampleDiscreteState(parts[i], parts[j]);
//                if (state == 1) {
//                    bindingsListPatch1[i] += 1;
//                    bindingsListPatch2[j] += 1;
//                } else if (state == 2) {
//                    bindingsListPatch1[i] += 1;
//                    bindingsListPatch1[j] += 1;
//                } else if (state == 3) {
//                    bindingsListPatch2[i] += 1;
//                    bindingsListPatch1[j] += 1;
//                } else if (state == 4) {
//                    bindingsListPatch2[i] += 1;
//                    bindingsListPatch2[j] += 1;
//                }
//            }
//        }
//        // Check conditions for bound patches are satisfied
//        for (int i = 0; i < parts.size(); i++) {
//            if (bindingsListPatch1[i] == 1) {
//                conditionBoundPatch1[i] = true;
//            }
//            if (bindingsListPatch2[i] == 1) {
//                conditionBoundPatch2[i] = true;
//            }
//        }
//        // Check if pentamer has formed
//        if (conditionBoundPatch1 == referenceCondition and conditionBoundPatch2 == referenceCondition) {
//            pentamerFormed = true;
//        }
//        return pentamerFormed;
//    }



}
