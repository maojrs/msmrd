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
        double relativeDistanceCutOff;
        int numParticleTypes;
    public:
        eventManager eventMgr = eventManager();
        msmrdMSM markovModel;
        spherePartition *positionPart;
        fullPartition *positionOrientationPart;

        /**
        * @param relativeDistanceCutOff radial cut off that define distance at which MSM/RD coupling is activated. This
        * parameter has to be consistent with the relative distance cut off used to obtain the discrete trajectories
        * (see patchyDimer trajectory for example).
        * @param numParticleTypes define the number of particle types in the unbound states in the
        * simulation (usually 1 or 2).
        * @param eventManager class to manage order of events (reactions/transitions).
        * @param markovModel pointer to class msmrdMarkovModel, which is the markovModel class specialized for
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
        msmrdIntegrator(double dt, long seed, std::string particlesbodytype, int numParticleTypes,
                        double relativeDistanceCutOff, std::vector<templateMSM> MSMlist, msmrdMSM markovModel);

        msmrdIntegrator(double dt, long seed, std::string particlesbodytype, int numParticleTypes,
                        double relativeDistanceCutOff, templateMSM MSMlist, msmrdMSM markovModel);

        // Redefine integrate function
        void integrate(std::vector<particleMS> &parts);

        void setDiscretization();
        void setDiscretization(spherePartition *thisPositionPartition);
        void setDiscretization(fullPartition *thisFullPartition);

        void computeTransitions2BoundStates(std::vector<particleMS> &parts);

        void computeTransitions2UnboundStates(std::vector<particleMS> &parts);

        void removeUnrealizedEvents(std::vector<particleMS> &parts);

        void transition2BoundState(std::vector<particleMS> &parts, int iIndex, int jIndex, int endState);

        void transition2UnboundState(std::vector<particleMS> &parts, int iIndex, int jIndex, int endState);

        float getRateFromKey(std::string);


    };

    /* Templated declarations (need to be in header). */

    /**
     * Constructors for MSM/RD integration class, can take discrete time MSM (msm) or a continous-time MSM.
     */
    template <typename templateMSM>
    msmrdIntegrator<templateMSM>::msmrdIntegrator(double dt, long seed, std::string particlesbodytype,
                                            int numParticleTypes, double relativeDistanceCutOff,
                                            std::vector<templateMSM> MSMlist, msmrdMSM markovModel) :
            overdampedLangevinMarkovSwitch<templateMSM>(MSMlist, dt, seed, particlesbodytype), markovModel(markovModel),
            numParticleTypes(numParticleTypes), relativeDistanceCutOff(relativeDistanceCutOff) {

        setDiscretization();

        if (MSMlist.size() != numParticleTypes and MSMlist.size() != 1) {
            std::__throw_range_error("Number of MSMs provided in MSMlist should match number of unbound particle"
                                     "types in the simulation. If same all particles share same MSM, is enough to "
                                     "provide one).");
        }
    };

    template <typename templateMSM>
    msmrdIntegrator<templateMSM>::msmrdIntegrator(double dt, long seed, std::string particlesbodytype,
                                            int numParticleTypes, double relativeDistanceCutOff,
                                                  templateMSM MSMlist, msmrdMSM markovModel) :
            overdampedLangevinMarkovSwitch<templateMSM>(MSMlist, dt, seed, particlesbodytype), markovModel(markovModel),
            numParticleTypes(numParticleTypes), relativeDistanceCutOff(relativeDistanceCutOff) {

        setDiscretization();

    };


    /* Sets default discretization (partition). Sets both positionPart and positionOrientationPart in case either of
     * the two is required. Note only one will actually be used by code. */
    template <typename templateMSM>
    void msmrdIntegrator<templateMSM>::setDiscretization() {
        int numSphericalSectionsPos = 7;
        int numRadialSectionsQuat = 5;
        int numSphericalSectionsQuat = 7;
        positionPart = new spherePartition(numSphericalSectionsPos);
        positionOrientationPart = new fullPartition(relativeDistanceCutOff, numSphericalSectionsPos,
                                                    numRadialSectionsQuat, numSphericalSectionsQuat);
    }

    // Sets pointer to discretization chosen (no rotation discretization).
    template <typename templateMSM>
    void msmrdIntegrator<templateMSM>::setDiscretization(spherePartition *thisPositionPartition) {
        delete positionPart;
        positionPart = thisPositionPartition;
    }

    // Sets pointer to discretization chosen (discretization taking into account rotation).
    template <typename templateMSM>
    void msmrdIntegrator<templateMSM>::setDiscretization(fullPartition *thisFullPartition) {
        delete positionOrientationPart;
        positionOrientationPart = thisFullPartition;
    }

    template <typename templateMSM>
    float msmrdIntegrator<templateMSM>::getRateFromKey(std::string key) {
        return markovModel.getRate(key);
    };



}
