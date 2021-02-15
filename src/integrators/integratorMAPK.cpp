//
// Created by maojrs on 2/15/21.
//

#include "integrators/integratorMAPK.hpp"

namespace msmrd {

    integratorMAPK::integratorMAPK(double dt, long seed, std::string particlesbodytype,
                                   int reactivationRateK, int reactivationRateP, std::vector<int> mapkIndex,
                                   std::vector<int> kinaseIndex, std::vector<int> phosIndex) :
                                   overdampedLangevin(dt, seed, particlesbodytype),
                                   reactivationRateK(reactivationRateK), reactivationRateP(reactivationRateP),
                                   mapkIndex(mapkIndex), kinaseIndex(kinaseIndex), phosIndex(phosIndex) {}

    // Integrate one particle from the particle list main routine (visible only inside the class)
    void integratorMAPK::integrateOne(int partIndex, std::vector<particle> &parts, double timestep) {
        vec3<double> force;
        vec3<double> torque;
        force = 1.0*forceField[partIndex];
        torque = 1.0*torqueField[partIndex];
        translate(parts[partIndex], force, timestep);
        if (rotation) {
            rotate(parts[partIndex], torque, timestep);
        }
        assignState(partIndex, parts);
    }

    void integratorMAPK::assignState(int partIndex, std::vector<particle> &parts) {
        // If current particle corresponds to an MAPK molecule, check binding
        double reactivationRate = 0;
        std::array<int, 2> MAPKbindings;
        int stateMAPK;
        if (std::find(mapkIndex.begin(), mapkIndex.end(), partIndex) != mapkIndex.end()) {
            MAPKbindings = checkMAPKbindings(partIndex, parts);
            // Assign bound states to kinase and phosphate molecules and time to reactivation
            for (auto &bindingIndex : MAPKbindings) {
                if (bindingIndex > -1) {
                    parts[bindingIndex].state = 1;
                    if (std::find(kinaseIndex.begin(), kinaseIndex.end(), bindingIndex) != kinaseIndex.end()){
                        reactivationRate = reactivationRateK;
                    } else {
                        reactivationRate = reactivationRateP;
                    }
                    auto r1 = randg.uniformRange(0, 1);
                    parts[bindingIndex].timeCounter = std::log(1.0 / r1) / reactivationRate;
                }
            }
            // Assign state to MAPK molecules
            if (MAPKbindings[0] == -1 && MAPKbindings[0] == -1){
                parts[partIndex].state = 0;
            } else if (MAPKbindings[0] > -1 && MAPKbindings[0] == -1) {
                parts[partIndex].state = 1;
            } else if (MAPKbindings[0] > -1 && MAPKbindings[0] == -1) {
                parts[partIndex].state = 2;
            } else{
                parts[partIndex].state = 3;
            }
        }
    }

    /* Returns tuple of integers corresponding to indexes of particle to which it is bound. If first
    * value is -1, particle is not bound on that site. This function only checks for bindings in the
    * MAPK molecules */
    std::array<int, 2> integratorMAPK::checkMAPKbindings(int partIndex, std::vector<particle> &parts) {
        int bindingIndexAt1 = -1;
        int bindingIndexAt2 = -1;
        std::array<int, 2> bindings;
        // NEED TO WRITE ALGORITHM TO FIND BINDINGS
        bindings = {bindingIndexAt1, bindingIndexAt2};
        return bindings;
    }


}
