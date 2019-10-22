//
// Created by maojrs on 10/7/19.
//

#pragma once

#include "integrators/msmrdIntegrator.hpp"

namespace msmrd {

    using msm = msmrd::discreteTimeMarkovStateModel;
    using ctmsm = msmrd::continuousTimeMarkovStateModel;
    using msmrdMSM = msmrd::msmrdMarkovModel;
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

        void computeTransitionsFromTransitionStates(std::vector<particleMS> &parts) override;

        void computeTransitionsFromBoundStates(std::vector<particleMS> &parts) override;

        void transition2BoundState(std::vector<particleMS> &parts, int iIndex, int jIndex, int endState) override;

        //void transition2UnboundState(std::vector<particleMS> &parts, int iIndex, int jIndex, int endState) override;

        //void transitionBetweenBoundStates(std::vector<particleMS> &parts, int iIndex,
        //                                  int jIndex, int endState) override;

        //void removeUnrealizedEvents(std::vector<particleMS> &parts) override;

    protected:
        /* Functions exclusive to multiparticle MSM/RD*/

        int addComplex(std::vector<particleMS> &parts, int iIndex, int jIndex, int endState);

        void updateParticleComplexesVector(std::vector<particleMS> &parts);

        void setCompoundPositionOrientation(std::vector<particleMS> &parts, int iIndex, int jIndex,
                                            int mainComplexSize);


    };

};

