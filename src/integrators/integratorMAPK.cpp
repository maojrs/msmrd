//
// Created by maojrs on 2/15/21.
//

#include "integrators/integratorMAPK.hpp"

namespace msmrd {

    integratorMAPK::integratorMAPK(double dt, long seed, std::string particlesbodytype, double anglePatches,
                                   double reactivationRateK, double reactivationRateP, std::vector<int> mapkIndex,
                                   std::vector<int> kinaseIndex, std::vector<int> phosIndex) :
                                   overdampedLangevin(dt, seed, particlesbodytype), anglePatches(anglePatches),
                                   reactivationRateK(reactivationRateK), reactivationRateP(reactivationRateP),
                                   mapkIndex(mapkIndex), kinaseIndex(kinaseIndex), phosIndex(phosIndex) {
        // Define patches locations in terms of the anglePatches variable
        MAPKpatch1 = {std::cos(anglePatches/2.0),std::sin(anglePatches/2.0), 0.};
        MAPKpatch2 = {std::cos(-anglePatches/2.0),std::sin(-anglePatches/2.0), 0.};
        ligandPatch = {std::cos(-anglePatches/2.0),std::sin(-anglePatches/2.0), 0.};
//        // Define Markov models for reactivation of kinase and phosphotase
//        int msmidKinase = 1;
//        int msmidPhos = 2
//        long seed = -1;
//        std::vector<std::vector<double>> tmatrixKinase ={ {-reactivationRateK, reactivationRateK},
//                                                    {reactivationRateK, -reactivationRateK} };
//        std::vector<std::vector<double>> tmatrixPhos ={ {-reactivationRateP, reactivationRateP},
//                                                          {reactivationRateP, -reactivationRateP} };
//        ctmsmKinase = std::make_unique<ctmsm>(msmidKinase, tmatrixKinase, seed);
//        ctmsmPhos = std::make_unique<ctmsm>(msmidPhos, tmatrixPhos, seed);

    }

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
        reactivationKorP(partIndex, parts, timestep);
        assignState(partIndex, parts);
    }

    void integratorMAPK::assignState(int partIndex, std::vector<particle> &parts) {
        double reactivationRate = 0;
        int stateMAPK;
        std::array<bool,2> MAPKphosphorilated = {false, false}; // site1 and site2 respectively
        // If current particle corresponds to an MAPK molecule, check binding and assign states
        if (std::find(mapkIndex.begin(), mapkIndex.end(), partIndex) != mapkIndex.end()) {
            // Check current phosphorilated state of MAPK
            if (parts[partIndex].state == 1 or parts[partIndex].state == 3){
                MAPKphosphorilated[0] = true;
            }
            if (parts[partIndex].state == 2 or parts[partIndex].state == 3){
                MAPKphosphorilated[1] = true;
            }
            // Check for current bindings
            auto MAPKbindingsAndTypes = checkMAPKbindings(partIndex, parts);
            // Assign bound states to kinase and phosphate molecules and time to reactivation
            for (int i = 0; i < MAPKbindingsAndTypes.size(); i++) {
                auto bindingIndex = std::get<0>(MAPKbindingsAndTypes[i]);
                auto bindingType = std::get<1>(MAPKbindingsAndTypes[i]);
                if (bindingIndex > -1) {
                    parts[bindingIndex].state = 1; // set ligand particle to inactive state
                    if (bindingType == 1){
                        if (not MAPKphosphorilated[i]) {
                            reactivationRate = 1.0 * reactivationRateK;
                            MAPKphosphorilated[i] = true;
                        }
                    } else { //bindingType 2
                        if (MAPKphosphorilated[i]) {
                            reactivationRate = 1.0 * reactivationRateP;
                            MAPKphosphorilated[i] = false;
                        }
                    }
                    auto r1 = randg.uniformRange(0, 1);
                    parts[bindingIndex].timeCounter = std::log(1.0 / r1) / reactivationRate;
                }
            }
            auto MAPKbindings = std::array<int,2>{std::get<0>(MAPKbindingsAndTypes[0]),
                                                  std::get<0>(MAPKbindingsAndTypes[1])};
            auto MAPKtypes = std::array<int,2>{std::get<1>(MAPKbindingsAndTypes[0]),
                                               std::get<1>(MAPKbindingsAndTypes[1])};
            // Assign state to MAPK molecules
            if (not MAPKphosphorilated[0] and not MAPKphosphorilated[1]){
                parts[partIndex].state = 0;
            } else if (MAPKphosphorilated[0] && not MAPKphosphorilated[1]) {
                parts[partIndex].state = 1;
            } else if (not MAPKphosphorilated[0] && MAPKphosphorilated[1]) {
                parts[partIndex].state = 2;
            } else{
                parts[partIndex].state = 3;
            }
        }
    }

    /* Checks if the current molecule, in case of being a kinase or phosphotase, needs to be reactivated. If so,
     * it reacivates the molecule. It also updates the time counter for inactive particles */
    void integratorMAPK::reactivationKorP(int partIndex, std::vector<particle> &parts, double timestep){
        if (parts[partIndex].type >0 and parts[partIndex].state == 1) {
            parts[partIndex].timeCounter -= timestep;
            if (parts[partIndex].timeCounter <= 0){
                // Reactivate particle
                parts[partIndex].state = 0;
            }
        }
    }

    /* Returns array of tuples. The first tuple corresponds to index and type bound to site 1. If index value
     * is -1, particle is not bound on that site. Te second tuple the same for site 2 The type cn be 1 or 2
     * depending if the kinase or phosphotase are bound. Note we only checks for binding in the MAPK molecules. */
    std::array<std::tuple<int,int>, 2> integratorMAPK::checkMAPKbindings(int partIndex,
            std::vector<particle> &parts) {
        int bindingIndexAt1 = -1;
        int bindingIndexAt2 = -1;
        int bindingTypeAt1 = -1;
        int bindingTypeAt2 = -1;
        bool boundAt1 = false;
        bool boundAt2 = false;
        auto particleMAPK = &parts[partIndex];
        std::array<int, 2> bindings;
        std::array<int, 2> bindingTypes;
        vec3<double> relPosition;
        // Finds bindings following assumptions of patch locations in header file
        for (int i = 0; i < parts.size(); i++) {
            // If kinase or phosphotase
            if (parts[i].type > 0) {
                relPosition = calculateRelativePosition(particleMAPK->position, parts[i].position);
                relPosition = msmrdtools::rotateVec(relPosition, particleMAPK->orientation.conj());
                if ((MAPKpatch1 - relPosition).norm() <= toleranceBinding and not boundAt1) {
                    // Calculate vectors pointing to binding patches
                    auto rotatedMAPKpatch1 = msmrdtools::rotateVec(MAPKpatch1, particleMAPK->orientation);
                    auto rotatedLigandPatch = msmrdtools::rotateVec(ligandPatch, parts[i].orientation);
                    if ((rotatedMAPKpatch1 + rotatedLigandPatch).norm() <= toleranceBinding){
                        //BINDING
                        bindingIndexAt1 = i;
                        bindingTypeAt1 = parts[i].type;
                        boundAt1 = true;
                    }
                }
                if ((MAPKpatch2 - relPosition).norm() <= toleranceBinding and not boundAt2) {
                    // Calculate vectors pointing to binding patches
                    auto rotatedMAPKpatch2 = msmrdtools::rotateVec(MAPKpatch2, particleMAPK->orientation);
                    auto rotatedLigandPatch = msmrdtools::rotateVec(ligandPatch, parts[i].orientation);
                    if ((rotatedMAPKpatch2 + rotatedLigandPatch).norm() <= toleranceBinding){
                        //BINDING
                        bindingIndexAt2 = i;
                        bindingTypeAt2 = parts[i].type;
                        boundAt2 = true;
                    }
                }
            }
        }
        auto bindingAt1 = std::make_tuple(bindingIndexAt1, bindingTypeAt1);
        auto bindingAt2 = std::make_tuple(bindingIndexAt2, bindingTypeAt2);
        return {bindingAt1, bindingAt2};
    }


}
