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
     * Class for multi-particle msmrd integration based on patchy particles.
     */
    template <typename templateMSM>
    class msmrdMultiParticleIntegrator : public msmrdIntegratorDiscrete<ctmsm> {
    public:
        using msmrdIntegratorDiscrete<ctmsm>::msmrdIntegratorDiscrete;

        void transition2BoundState(std::vector<particleMS> &parts, int iIndex, int jIndex, int endState);
    };

};

