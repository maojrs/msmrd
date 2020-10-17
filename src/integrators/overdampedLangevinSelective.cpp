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
        // Integrate and save next positions/orientations in parts[i].next***
        for (int i = 0; i < parts.size(); i++) {
            setActivePatches(parts);
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

    /* Check if pentamer was formed
     */
    bool overdampedLangevinSelective::hasPentamerFormed(std::vector<particle> &parts) {
        auto pentamerFormed = false;
        std::vector<int> bindingsListPatch1(5, 0);
        std::vector<int> bindingsListPatch2(5, 0);
        // Count bindings in patches
        for (int i = 0; i < parts.size() - 1; i++) {
            for (int j = i + 1; j < parts.size(); j++) {
                // Only allow a new bound state if there is no previous bound state between
                auto state = patchyTraj->sampleDiscreteState(parts[i], parts[j]);
                if (state == 1) {
                    bindingsListPatch1[i] += 1;
                    bindingsListPatch2[j] += 1;
                } else if (state == 2) {
                    bindingsListPatch1[i] += 1;
                    bindingsListPatch1[j] += 1;
                } else if (state == 3) {
                    bindingsListPatch2[i] += 1;
                    bindingsListPatch1[j] += 1;
                } else if (state == 4) {
                    bindingsListPatch2[i] += 1;
                    bindingsListPatch2[j] += 1;
                }
            }
        }
        // Check conditions for bound patches are satisfied
        for (int i = 0; i < parts.size(); i++) {
            if (bindingsListPatch1[i] >= 1) {
                conditionBoundPatch1[i] = true;
            }
            if (bindingsListPatch2[i] >= 1) {
                conditionBoundPatch2[i] = true;
            }
        }
        // Check if pentamer has formed
        if (conditionBoundPatch1 == referenceCondition and conditionBoundPatch2 == referenceCondition) {
            pentamerFormed = true;
        }
        return pentamerFormed;
    }



}
