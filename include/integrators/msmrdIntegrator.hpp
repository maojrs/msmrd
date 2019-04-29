//
// Created by maojrs on 2/6/19.
//

#pragma once

#include "discretizations/positionOrientationPartition.hpp"
#include "eventManager.hpp"
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
        double relativeDistanceCutOff = 2.2;
        msmrdMSM &markovModel;
        fullPartition &positionOrientationPart;
        eventManager eventMgr = eventManager();
    public:
        /**
        * @param relativeDistanceCutOff radial cut off that define distance at which MSM/RD coupling is activated. This
        * parameter has to be consistent with the relative distance cut off used to obtain the discrete trajectories
        * (see patchyDimer trajectory for example).
        * @param markovModel pointer to class msmrdMarkovModel, which is the markovModel class specialized for
        * the MSM/RD scheme.
        * @param positionOrientationPart pointer to full partition of relative distance and relative orientation,
        * a.k.a positionOrientationPartition.
        * @param eventManager class to manage order of events (reactions/transitions).
        */
        msmrdIntegrator(double dt, long seed, std::string particlesbodytype, double relativeDistanceCutOff,
                        std::vector<templateMSM> MSMlist, msmrdMSM markovModel, fullPartition positionOrientationPart);

        // Redefine integrate function
        void integrate(std::vector<particleMS> &parts);

        void computeTransitions2BoundStates(std::vector<particleMS> &parts);

        void computeTransitions2UnboundStates(std::vector<particleMS> &parts);

    };




}
