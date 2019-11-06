//
// Created by maojrs on 11/5/19.
//

#pragma once

#include "integrators/msmrdIntegrator.hpp"

namespace msmrd {

    using msm = msmrd::discreteTimeMarkovStateModel;
    using ctmsm = msmrd::continuousTimeMarkovStateModel;
    using msmrdMSM = msmrd::msmrdMarkovModel;
    using fullPartition = msmrd::positionOrientationPartition;

    /**
     * Class for msmrd integration of patchy proteins. It is an specialization of the msmrdIntegrator class
     * that allows particle two to have two states/conformation. It is not yet implemented for a general case of
     * n and m states conformations of each particle. Currently it only m
     */
    class msmrdPatchyProtein : public msmrdIntegrator<ctmsm> {
    public:
        using msmrdIntegrator<ctmsm>::msmrdIntegrator;

    };

}
