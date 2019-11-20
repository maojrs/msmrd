//
// Created by maojrs on 8/19/19.
//


#pragma once

#include <utility>
#include "discretizations/positionOrientationPartition.hpp"
#include "integrators/overdampedLangevinMarkovSwitch.hpp"
#include "markovModels/msmrdMarkovModel.hpp"
#include "eventManager.hpp"
#include "tools.hpp"

namespace msmrd {

    using msm = msmrd::discreteTimeMarkovStateModel;
    using ctmsm = msmrd::continuousTimeMarkovStateModel;
    using msmrdMSM = msmrd::msmrdMarkovModel;
    using fullPartition = msmrd::positionOrientationPartition;

    /**
     * Base class for msmrd integration (coupling MSM and reaction-diffusion)
     * @tparam templateMSM template can be an msm or a ctmsm. This only referes to the
     * MSM used when the particles are unbound in MSM/RD. If they are bound, they use the only
     * MSM available at the moment for MSM/RD integration, which is based on a discrete-time MSM
     * obtained with PyEmma.
     */
    template <typename templateMSM>
    class msmrdIntegrator : public overdampedLangevinMarkovSwitch<templateMSM> {
    protected:
        std::array<double,2> radialBounds;
        int numParticleTypes;
        bool firstrun = true;
        bool recordEventLog = false;
    public:
        eventManager eventMgr = eventManager();
        msmrdMSM markovModel;
        //spherePartition *positionPart;
        //fullPartition *positionOrientationPart;
        std::shared_ptr<spherePartition> positionPart;
        std::shared_ptr<fullPartition> positionOrientationPart;

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
        * @param firstrun boolean variable to check if the integrator is ran for the first time in a simulation.
        * @param recordEventLog boolean to dump or not dump event list for every time step
        * into eventMgr.eventLog. Useful for debugging, but need to watch out memory if log not dumped fast enough.
        * @param eventManager class to manage order of events (reactions/transitions).
        * @param markovModel pointer to class msmrdMSMDiscrete, which is the markovModel class specialized for
        * the MSM/RD scheme. It controls the markov Model in the bound state and the msmrd coupling.
        * @param positionPart shared pointer to spherical partition of relative distance, only used when rotation
        * is disabled in the integrator. Needs to be shared so discretization can be set from python interface.
        * @param positionOrientationPart shared pointer to full partition of relative distance and relative orientation,
        * a.k.a positionOrientationPartition, used when rotation is enabled. Needs to be shared so discretization
        * can be set from python interface.
        * @param MSMlist (see overdampedLanegvinMarkovSwitch parent class) can be either a vector of msms
        * (not yet implemented) or of ctmsms. It corresponds to the MSMs of the unbound particles (conformation
        * switching). Its size should match the number of unbound particle types in the simulation, one MSM
        * per particle type. Size 1 is also acceptable, which assumes all aprticles share the same MSM.
        */
        msmrdIntegrator(double dt, long seed, std::string particlesbodytype, int numParticleTypes,
                        std::array<double,2> radialBound, std::vector<templateMSM> MSMlist,
                                msmrdMSM markovModel);

        msmrdIntegrator(double dt, long seed, std::string particlesbodytype, int numParticleTypes,
                        std::array<double,2> radialBound, templateMSM MSMlist, msmrdMSM markovModel);


        // Redefine integrate function
        void integrate(std::vector<particle> &parts) override;


        // Auxiliary functions for main MSM/RD functions (can be set to virtual if they need to be overriden)
        void integrateDiffusion(std::vector<particle> &parts, double dt);

        std::tuple<vec3<double>, quaternion<double>> getRelativePositionOrientation(int state);

        int setRandomUnboundState(std::vector<particle> &parts, int partIndex);

        void setDefaultDiscretization();

        void setDiscretization(std::shared_ptr<spherePartition> &thisSpherePartition);

        void setDiscretization(std::shared_ptr<fullPartition> &thisFullPartition);


        /* Other useful functions to output log files (setRecordEventLog sets recordEventLog variable. This
         * needs to be set to true in order to printEventLog function to print anything meaningful) */

        void setRecordEventLog(bool value) { recordEventLog = value; }

        void printEventLog(std::string filename);


        /* Main MSM/RD function. They are defined as virtual in case we want to override them in derived classes
         * to modify fucntionality */

        virtual int computeCurrentTransitionState(particle &part1, particle &part2);

        virtual void computeTransitionsFromTransitionStates(std::vector<particle> &parts);

        virtual void computeTransitionsFromBoundStates(std::vector<particle> &parts);

        virtual void transition2BoundState(std::vector<particle> &parts, int iIndex, int jIndex, int endState);

        virtual void transition2UnboundState(std::vector<particle> &parts, int iIndex, int jIndex, int endState);

        virtual void transitionBetweenBoundStates(std::vector<particle> &parts, int iIndex, int jIndex, int endState);

        virtual void transitionBetweenTransitionStates(int iIndex, int jIndex);

        virtual void removeUnrealizedEvents(std::vector<particle> &parts);

        virtual void applyEvents(std::vector<particle> &parts);




    };

    /* Templated declarations (need to be in header). */

    /**
     * Constructors for MSM/RD integration class, can take discrete time MSM (msm) or a continous-time MSM.
     */
    template <typename templateMSM>
    msmrdIntegrator<templateMSM>::msmrdIntegrator(double dt, long seed,
                               std::string particlesbodytype, int numParticleTypes, std::array<double,2> radialBounds,
                               std::vector<templateMSM> MSMlist, msmrdMSM markovModel) :
            overdampedLangevinMarkovSwitch<templateMSM>(MSMlist, dt, seed, particlesbodytype),
                    markovModel(std::move(markovModel)),
            numParticleTypes(numParticleTypes), radialBounds(radialBounds) {

        setDefaultDiscretization();

        if (MSMlist.size() != numParticleTypes and MSMlist.size() != 1) {
            throw std::invalid_argument("Number of MSMs provided in MSMlist should match number of unbound particle"
                                        "types in the simulation. If same all particles share same MSM, is enough to "
                                        "provide one).");
        }
    };

    template <typename templateMSM>
    msmrdIntegrator<templateMSM>::msmrdIntegrator(double dt, long seed,
                               std::string particlesbodytype, int numParticleTypes, std::array<double,2> radialBounds,
                               templateMSM MSMlist, msmrdMSM markovModel) :
            overdampedLangevinMarkovSwitch<templateMSM>(MSMlist, dt, seed, particlesbodytype),
                    markovModel(std::move(markovModel)),
            numParticleTypes(numParticleTypes), radialBounds(radialBounds) {

        setDefaultDiscretization();

    };


    /* Sets default discretization (partition). Sets both positionPart and positionOrientationPart in case either of
     * the two is required. Note only one will actually be used by code. */
    template <typename templateMSM>
    void msmrdIntegrator<templateMSM>::setDefaultDiscretization() {
        int numSphericalSectionsPos = 7;
        int numRadialSectionsQuat = 5;
        int numSphericalSectionsQuat = 7;
        double relativeDistanceCutOff = radialBounds[1];

        positionPart = std::make_shared<spherePartition>(spherePartition(numSphericalSectionsPos));
        positionOrientationPart = std::make_shared<fullPartition> (fullPartition(relativeDistanceCutOff,
                numSphericalSectionsPos, numRadialSectionsQuat, numSphericalSectionsQuat));
    }


    // Sets pointer to discretization chosen (discretization spherePartition, no rotation).
    template <typename templateMSM>
    void msmrdIntegrator<templateMSM>::setDiscretization(std::shared_ptr<spherePartition> &thisSpherePartition) {
        positionPart.reset();
        positionPart = thisSpherePartition;
    }

    // Sets pointer to discretization chosen (fullPartition discretization taking into account rotation).
    template <typename templateMSM>
    void msmrdIntegrator<templateMSM>::setDiscretization(std::shared_ptr<fullPartition> &thisFullPartition) {
        positionOrientationPart.reset();
        positionOrientationPart = thisFullPartition;
    }

    // Prints eventlog by invoking method from eventMgr into file filename.dat
    template <typename templateMSM>
    void msmrdIntegrator<templateMSM>::printEventLog(std::string filename) {
        if (not recordEventLog) {
            throw std::runtime_error("Need to use setRecordEventLog to set recordEventLog true in order to "
                                     "print the event log.");
        }
        eventMgr.printEventLog(filename);
    }

}
