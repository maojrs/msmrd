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
     * n and m states conformations of each particle.
     */
    class msmrdPatchyProtein : public msmrdIntegrator<ctmsm> {
    public:
        using msmrdIntegrator<ctmsm>::msmrdIntegrator;

        // Auxiliary functions

        int setNewUnboundState(int unboundState, std::vector<particle> &parts, int partIndex);


        // Main overridden functions

        int computeCurrentTransitionState(particle &part1, particle &part2) override;

        void transition2UnboundState(std::vector<particle> &parts, int iIndex, int jIndex, int endState) override;



    };

}
