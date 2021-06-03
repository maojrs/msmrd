//
// Created by maojrs on 5/31/21.
//

#pragma once

#include "integrators/msmrdMultiParticleIntegrator.hpp"

namespace msmrd {

    using ctmsm = msmrd::continuousTimeMarkovStateModel;
    using fullPartition = msmrd::positionOrientationPartition;

    /**
     * Class for multi-particle msmrd integration for the MAPK application. Uses base functionality from
     * msmrdMultiParticleIntegrator, but changes some methods to adapt to the MAPK application.
     */
    class msmrdMAPKintegrator : public msmrdMultiParticleIntegrator<ctmsm> {
    public:
        using msmrdMultiParticleIntegrator::msmrdMultiParticleIntegrator;
    };


}
