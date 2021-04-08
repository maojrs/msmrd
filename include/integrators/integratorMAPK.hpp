//
// Created by maojrs on 2/15/21.
//

#pragma once
#include "integrators/overdampedLangevin.hpp"
#include "particle.hpp"
#include "potentials/potentials.hpp"
#include "markovModels/continuousTimeMarkovModel.hpp"

using ctmsm = msmrd::continuousTimeMarkovStateModel;

namespace msmrd {
    /**
     * Integrator class for integratorMAPK pathway simulation. In addition to the usual overdamped Langevin dynamics
     * it checks for the state and if phosphorilation happened. Uses the particle type to know if it is an
     * integratorMAPK molecule (type 0), a Kinase(type 1) or a phosphatase (type 2).
     *
     * We further use states to determine the current status of each molecule:
     * integratorMAPK state 0: neither of the two binding sites has been phosphorylated
     * integratorMAPK state 1: only binding site 1 is phosphorylated.
     * integratorMAPK state 2: only binding site 2 is phosphorylated.
     * integratorMAPK state 3: both binding sites are phosphorylated.
     * Kinase state 0: active
     * Kinase state 1: inactive
     * Phosphatase state 0: active
     * Phosphatase state 1: inactive
     */
    class integratorMAPK : public overdampedLangevin {
    protected:
        double anglePatches;
        double reactivationRateK;
        double reactivationRateP;
        std::vector<int> mapkIndex;
        std::vector<int> kinaseIndex;
        std::vector<int> phosIndex;
        double toleranceBinding = 0.2; //0.12;
        double toleranceBindingOrientation = 0.2;
        vec3<double> MAPKpatch1;
        vec3<double> MAPKpatch2;
        vec3<double> ligandPatch;
//        std::unique_ptr<ctmsm> ctmsmKinase;
//        std::unique_ptr<ctmsm> ctmsmPhos;
        /**
         * @param anglePatches angle between binding sites (patches) os MAPK molecules. The patches are
         * assumed to be located at anglePatches/2 and -anglePacthes/2 for the MAPK molecules. And
         * at -anglePatches/2 for the one patch of kinases and phosphotases. (Both in the xy plane).
         * @param reactivationRateK reactivation rate for kinase after binding
         * @param reactivationRateP reactivation rate for phosphatase after binding
         * @param mapkIndex is the indexes in the particle list of particles that correspond to
         * MAPK molecules
         * @param kinaseIndex same as mapkIndex but for kinase molecules
         * @param phosIndex same as mapkIndex but for phosphatase molecules
         * @param toleranceBinding tolerance for evaluating if relative position corresponds
         * to binding configuration
         * @param toleranceBinding tolerance for evaluating if relative orientation corresponds
         * to binding configuration
         * @param MAPKpatch1/2 unit vector pointing to patch/binding site 1 or 2 of MAPK molecule in its initial
         * orientation (1,0,0,0).
         * @param ligandPatch same as above for the ligand (kinase or phosphotase) patch.
         */

        void integrateOne(int partIndex, std::vector<particle> &parts, double timestep) override;

        std::array<std::tuple<int,int>, 2> checkMAPKbindings(int partIndex, std::vector<particle> &parts);

        void assignState(int partIndex, std::vector<particle> &parts);

        void reactivationKorP(particle &part, double timestep);
    public:
        integratorMAPK(double dt, long seed, std::string particlesbodytype,
                double anglePatches, double reactivationRateK, double reactivationRateP, std::vector<int> mapKIndex,
                std::vector<int> kinaseIndex,std::vector<int> phosIndex);
    };

}
