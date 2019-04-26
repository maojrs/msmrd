//
// Created by maojrs on 2/6/19.
//

#pragma once
#include "discretizations/positionOrientationPartition.hpp"
#include "integrators/overdampedLangevinMarkovSwitch.hpp"
#include "markovModels/msmrdMarkovModel.hpp"
#include "tools.hpp"

namespace msmrd {
    using msm = msmrd::discreteTimeMarkovStateModel;
    using ctmsm = msmrd::continuousTimeMarkovStateModel;
    using msmrdMSM = msmrd::msmrdMarkovStateModel;
    using fullPartition = msmrd::positionOrientationPartition;

    /**
     * Base class for msmrd integration (coupling MSM and reaction-diffusion)
     * @tparam templateMSM template can be an msm or a ctmsm
     */
    template <typename templateMSM>
    class msmrdIntegrator : public overdampedLangevinMarkovSwitch<templateMSM> {
    private:
        msmrdMSM &markovModel;
        fullPartition &positionOrientationPart;
        std::vector<particleMS> boundParticles{};

    public:
        /**
        * @param markovModel pointer to class msmrdMarkovModel, which is the markovModel class specialized for
        * the MSM/RD scheme.
        * @param positionOrientationPart pointer to full partition of relative distance and relative orientation,
        * a.k.a positionOrientationPartition.
        * @param boundParticles vector of particles, where each particle represents the bound state between
         * two particles from the main particle list being integrated.
        */
        msmrdIntegrator(double dt, long seed, std::string particlesbodytype, std::vector<templateMSM> MSMlist,
                msmrdMSM markovModel, fullPartition positionOrientationPart);

        // Redefine integrate function
        void integrate(std::vector<particleMS> &parts);
    };




}
