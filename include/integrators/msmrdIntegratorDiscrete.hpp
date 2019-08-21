//
// Created by maojrs on 8/19/19.
//


#pragma once

#include <utility>
#include "discretizations/positionOrientationPartition.hpp"
#include "eventManager.hpp"
#include "integrators/overdampedLangevinMarkovSwitch.hpp"
#include "markovModels/msmrdMarkovModelDiscrete.hpp"
#include "tools.hpp"

namespace msmrd {

    using msm = msmrd::discreteTimeMarkovStateModel;
    using ctmsm = msmrd::continuousTimeMarkovStateModel;
    using msmrdMSMDiscrete = msmrd::msmrdMarkovModelDiscrete;
    using fullPartition = msmrd::positionOrientationPartition;

    /**
     * Base class for msmrd integration (coupling MSM and reaction-diffusion)
     * @tparam templateMSM template can be an msm or a ctmsm
     */
    template <typename templateMSM>
    class msmrdIntegratorDiscrete : public overdampedLangevinMarkovSwitch<templateMSM> {
    private:
        std::array<double,2> radialBounds;
        int numParticleTypes;
        bool firstrun = true;
    public:
        eventManager eventMgr = eventManager();
        msmrdMSMDiscrete markovModel;
        spherePartition *positionPart;
        fullPartition *positionOrientationPart;

        /**
        * @param radialBounds radial lower bound and upper bound of transition region in MSM/RD scheme. It MUST
        * match the parameters used in discretization to obtain the MSM and the rateDictionary (see patchyDimer
        * trajectory for example). If relative distance is smaller than lower radial bound (radialBound[0]), we
        * are on the bound states region where. Between the bounds radialBound[0] and radialBound[1] the discretization
        * defines the transition states region. If the relative distance is larger than the upper radial bound
        * radialBound[1], then MSM/RD is completely desactivated and the dynamics are independent. The upper radial
        * bound radialBound[1], corresponds to the relativeDistanceCutOff in other sections of the code.
        * the discretization used coreMSM
        * @param numParticleTypes define the number of particle types in the unbound states in the
        * simulation (usually 1 or 2).
        * @param eventManager class to manage order of events (reactions/transitions).
        * @param markovModel pointer to class msmrdMSMDiscrete, which is the markovModel class specialized for
        * the MSM/RD scheme. It controls the markov Model in the bound state and the msmrd coupling.
        * @param positionPart pointer to spherical partition of relative distance, only used when rotation
        * is disabled in the integrator
        * @param positionOrientationPart pointer to full partition of relative distance and relative orientation,
        * a.k.a positionOrientationPartition, used when rotation is enabled.
        * @param MSMlist (see overdampedLanegvinMarkovSwitch parent class) can be either a vector of msms
        * (not yet implemented) or of ctmsms. It corresponds to the MSMs of the unbound particles (conformation
        * switching). Its size should match the number of unbound particle types in the simulation, one MSM
        * per particle type. Size 1 is also acceptable, which assumes all aprticles share the same MSM.
        */
        msmrdIntegratorDiscrete(double dt, long seed, std::string particlesbodytype, int numParticleTypes,
                        std::array<double,2> radialBound, std::vector<templateMSM> MSMlist,
                                msmrdMSMDiscrete markovModel);

        msmrdIntegratorDiscrete(double dt, long seed, std::string particlesbodytype, int numParticleTypes,
                        std::array<double,2> radialBound, templateMSM MSMlist, msmrdMSMDiscrete markovModel);

        // Redefine integrate function
        void integrate(std::vector<particleMS> &parts);

        void integrateDiffusion(std::vector<particleMS> &parts, double dt);

        void setDefaultDiscretization();

        void setDiscretization(fullPartition *thisFullPartition);

        int computeCurrentTransitionState(particleMS &part1, particleMS &part2);

        void computeTransitionsFromTransitionStates(std::vector<particleMS> &parts);

        void computeTransitionsFromBoundStates(std::vector<particleMS> &parts);

        void transition2BoundState(std::vector<particleMS> &parts, int iIndex, int jIndex, int endState);

        void transition2UnboundState(std::vector<particleMS> &parts, int iIndex, int jIndex, int endState);

        void transitionBetweenBoundStates(std::vector<particleMS> &parts, int iIndex, int jIndex, int endState);

        void transitionBetweenTransitionStates(int iIndex, int jIndex);

        void removeUnrealizedEvents(std::vector<particleMS> &parts);

        void applyEvents(std::vector<particleMS> &parts);


    };

    /* Templated declarations (need to be in header). */

    /**
     * Constructors for MSM/RD integration class, can take discrete time MSM (msm) or a continous-time MSM.
     */
    template <typename templateMSM>
    msmrdIntegratorDiscrete<templateMSM>::msmrdIntegratorDiscrete(double dt, long seed,
                               std::string particlesbodytype, int numParticleTypes, std::array<double,2> radialBounds,
                               std::vector<templateMSM> MSMlist, msmrdMSMDiscrete markovModel) :
            overdampedLangevinMarkovSwitch<templateMSM>(MSMlist, dt, seed, particlesbodytype),
                    markovModel(std::move(markovModel)),
            numParticleTypes(numParticleTypes), radialBounds(radialBounds) {

        setDefaultDiscretization();

        if (MSMlist.size() != numParticleTypes and MSMlist.size() != 1) {
            std::__throw_range_error("Number of MSMs provided in MSMlist should match number of unbound particle"
                                     "types in the simulation. If same all particles share same MSM, is enough to "
                                     "provide one).");
        }
    };

    template <typename templateMSM>
    msmrdIntegratorDiscrete<templateMSM>::msmrdIntegratorDiscrete(double dt, long seed,
                               std::string particlesbodytype, int numParticleTypes, std::array<double,2> radialBounds,
                               templateMSM MSMlist, msmrdMSMDiscrete markovModel) :
            overdampedLangevinMarkovSwitch<templateMSM>(MSMlist, dt, seed, particlesbodytype),
                    markovModel(std::move(markovModel)),
            numParticleTypes(numParticleTypes), radialBounds(radialBounds) {

        setDefaultDiscretization();

    };


    /* Sets default discretization (partition). Sets both positionPart and positionOrientationPart in case either of
     * the two is required. Note only one will actually be used by code. */
    template <typename templateMSM>
    void msmrdIntegratorDiscrete<templateMSM>::setDefaultDiscretization() {
        int numSphericalSectionsPos = 7;
        int numRadialSectionsQuat = 5;
        int numSphericalSectionsQuat = 7;
        double relativeDistanceCutOff = radialBounds[1];
        positionPart = new spherePartition(numSphericalSectionsPos);
        positionOrientationPart = new fullPartition(relativeDistanceCutOff, numSphericalSectionsPos,
                                                    numRadialSectionsQuat, numSphericalSectionsQuat);
    }

    // Sets pointer to discretization chosen (discretization taking into account rotation).
    template <typename templateMSM>
    void msmrdIntegratorDiscrete<templateMSM>::setDiscretization(fullPartition *thisFullPartition) {
        delete positionOrientationPart;
        positionOrientationPart = thisFullPartition;
    }


}
