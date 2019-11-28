//
// Created by maojrs on 10/7/19.
//

#pragma once

#include "integrators/msmrdIntegrator.hpp"
#include "particleCompound.hpp"

namespace msmrd {

    using msm = msmrd::discreteTimeMarkovStateModel;
    using ctmsm = msmrd::continuousTimeMarkovStateModel;
    using fullPartition = msmrd::positionOrientationPartition;

    /**
     * Class for multi-particle msmrd integration based on patchy particles. Uses base functionality from
     * msmrdIntegratorDiscrete, but extends its methods for multiparticle MSM/RD integration.
     */
    class msmrdMultiParticleIntegrator : public msmrdIntegrator<ctmsm> {
    public:
        std::vector<particleCompound> particleCompounds;
        /**
         * @param particleCompounds: vector containing particle compounds (two or more particles bound
         * together) to track them and diffuse them.
         */

        using msmrdIntegrator<ctmsm>::msmrdIntegrator;

        void computeTransitionsFromTransitionStates(std::vector<particle> &parts) override;

        void computeTransitionsFromBoundStates(std::vector<particle> &parts) override;

        void transition2BoundState(std::vector<particle> &parts, int iIndex, int jIndex, int endState) override;

        void transition2UnboundState(std::vector<particle> &parts, int iIndex, int jIndex, int endState) override;

        //void transitionBetweenBoundStates(std::vector<particle> &parts, int iIndex,
        //                                  int jIndex, int endState) override;

        //void removeUnrealizedEvents(std::vector<particle> &parts) override;

    protected:

        /* Functions exclusive to multiparticle MSM/RD*/

        int addCompound(std::vector<particle> &parts, int iIndex, int jIndex, int endState);

        void setCompoundPositionOrientation(std::vector<particle> &parts, int iIndex, int jIndex,
                                            int mainComplexSize);

//        std::tuple<bool, bool> checkUnbindingCompounds(std::vector<particle> &parts, int iIndex, int jIndex);

        void updateParticleComplexesVector(std::vector<particle> &parts);

    };

};

