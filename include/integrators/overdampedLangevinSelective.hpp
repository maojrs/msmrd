//
// Created by maojrs on 10/15/20.
//
#pragma once
#include "integrators/overdampedLangevin.hpp"
#include "trajectories/discrete/discreteTrajectory.hpp"
#include "trajectories/discrete/patchyDimerTrajectory.hpp"



namespace msmrd {
    /**
     * Over-damped Langevin with patch selectivity integrator declaration. This
     * integrator is specialized for patchy particles, and it checks for active patches
     * to avoid for three particle bindings, i.e. one patch can only bind with one and only
     * one patch of another particle.
     *
     * CAREFUL: requires setting particles IDs corresponding to index of particle list and
     * setting particles activePatchList for every particle (see pentamerFPT script for example).
     * THIS is a specialized class for a specific application. For normal integration
     * use overdampedLangevin
     */
    class overdampedLangevinSelective : public overdampedLangevin {
    public:
        std::shared_ptr<discreteTrajectory<4>> patchyTraj;
        /**
         * @param discreteTraj pointer to patchy trajectory. The 4 indicate the number of bound states.
         * Can be extened to template to avoid hardcoding. As the functionality of this integrator is
         * heavily specialized, we live the number of bound states hard coded.
         */

        overdampedLangevinSelective(double dt, long seed, std::string particlesbodytype);

        void integrate(std::vector<particle> &parts) override;

        bool hasPentamerFormed(std::vector<particle> &parts);

        protected:

        void setActivePatches(std::vector<particle> &parts);
        };

}
