//
// Created by maojrs on 10/7/19.
//

#pragma once

#include "integrators/msmrdIntegratorDiscrete.hpp"

namespace msmrd {

    using msm = msmrd::discreteTimeMarkovStateModel;
    using ctmsm = msmrd::continuousTimeMarkovStateModel;
    using msmrdMSMDiscrete = msmrd::msmrdMarkovModelDiscrete;
    using fullPartition = msmrd::positionOrientationPartition;

    /**
     * Class for multi-particle msmrd integration based on patchy particles. Uses base functionality from
     * msmrdIntegratorDiscrete, but extends its methods for multiparticle MSM/RD integration.
     */
    class msmrdMultiParticleIntegrator : public msmrdIntegratorDiscrete<ctmsm> {
    public:
        std::vector<particleComplex> particleComplexes;
        /**
         * @param particleComplexes: vector containing particle complexex to track particles that are bound together.
         */

        using msmrdIntegratorDiscrete<ctmsm>::msmrdIntegratorDiscrete;

        void computeTransitionsFromTransitionStates(std::vector<particleMS> &parts) override;

        void computeTransitionsFromBoundStates(std::vector<particleMS> &parts) override;

        void transition2BoundState(std::vector<particleMS> &parts, int iIndex, int jIndex, int endState) override;

        //void transition2UnboundState(std::vector<particleMS> &parts, int iIndex, int jIndex, int endState) override;

        //void transitionBetweenBoundStates(std::vector<particleMS> &parts, int iIndex,
        //                                  int jIndex, int endState) override;

        //void removeUnrealizedEvents(std::vector<particleMS> &parts) override;

        /* Functions exclusive to multiparticle MSM/RD*/

        void addComplex(std::vector<particleMS> &parts, int iIndex, int jIndex, int endState);

        void updateParticleComplexesVector(std::vector<particleMS> &parts);


    };

};

