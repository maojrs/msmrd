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
     * Class for msmrd integration of patchy proteins. The normal MSM/RD scheme for patchyProtein trajectories
     * can be integrated with the msmrdIntegrator function directly. However using this integrator might
     * be more convenient since it redfines the default discretization.
     */
    class msmrdPatchyProtein : public msmrdIntegrator<ctmsm> {
    protected:

        void setPatchyProteinDiscretization();

    public:

        msmrdPatchyProtein(double dt, long seed, std::string particlesbodytype, int numParticleTypes,
                        std::array<double,2> radialBound, std::vector<ctmsm> MSMlist,
                        msmrdMSM markovModel);

        msmrdPatchyProtein(double dt, long seed, std::string particlesbodytype, int numParticleTypes,
                        std::array<double,2> radialBound, ctmsm MSMlist, msmrdMSM markovModel);
    };


    /* Specializations of this class are also included, specifically for the MSM/RD scheme derived from trajectories
     * based on patchyProteinTrajectory2. It is an specialization of ehe msmrdPatchyProtein class that takes into
     * account particle two having two states/conformation in its the discretization. */
    class msmrdPatchyProtein2 : public msmrdPatchyProtein {
    public:

        using msmrdPatchyProtein::msmrdPatchyProtein;

        // Auxiliary functions

        int setNewUnboundState(int unboundState, std::vector<particle> &parts, int partIndex);


        // Main overridden functions

        int computeCurrentTransitionState(particle &part1, particle &part2) override;

        void transition2UnboundState(std::vector<particle> &parts, int iIndex, int jIndex, int endState) override;


    };

}
