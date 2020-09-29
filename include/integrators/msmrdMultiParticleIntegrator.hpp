//
// Created by maojrs on 10/7/19.
//

#pragma once

#include "integrators/msmrdIntegrator.hpp"
#include "trajectories/discrete/discreteTrajectory.hpp"
#include "trajectories/discrete/patchyDimerTrajectory.hpp"
#include "particleCompound.hpp"

namespace msmrd {

    using msm = msmrd::discreteTimeMarkovStateModel;
    using ctmsm = msmrd::continuousTimeMarkovStateModel;
    using fullPartition = msmrd::positionOrientationPartition;

    /**
     * Class for multi-particle msmrd integration based on patchy particles. Uses base functionality from
     * msmrdIntegrator, but extends its methods for multiparticle MSM/RD integration.
     * @tparam templateMSM template can be an msm or a ctmsm. This only referes to the
     * MSM used when the particles are unbound in MSM/RD. If they are bound, they use the only
     * MSM available at the moment for MSM/RD integration, which is based on a discrete-time MSM
     * obtained with PyEmma. At the time, there is only an implementation using a ctmsm.
     */
    template <typename templateMSM>
    class msmrdMultiParticleIntegrator : public msmrdIntegrator<templateMSM> {
    public:
        std::vector<particleCompound> particleCompounds;
        std::shared_ptr<discreteTrajectory<4>> discreteTrajClass;
        const vec3<double> pentamerCenter = {1.0/(2.0 * std::sin(M_PI/5.0)), 0.0, 0.0};
        /**
         * @param particleCompounds: vector containing particle compounds (two or more particles bound
         * together) to track them and diffuse them.
         * @param discreteTrajClass: pointer to trajectory class used for the discretization. We need this to extract
         * the bound states and their relative positions and orientations. The specific discrete trajectory
         * class will be fixed in the constructor, so this needs to be modified for different MSM/RD applications.
         * @param pentamerCenter location of pentamer center with respect to the main particle in a compound. Uniquely
         * for trimer and pentamer implementations of multi-particle MSM/RD.
         */

        msmrdMultiParticleIntegrator(double dt, long seed, std::string particlesbodytype, int numParticleTypes,
                        std::array<double,2> radialBounds, std::vector<templateMSM> MSMlist,
                        msmrdMarkovModel msmrdMSM);

        msmrdMultiParticleIntegrator(double dt, long seed, std::string particlesbodytype, int numParticleTypes,
                        std::array<double,2> radialBounds, templateMSM MSMlist, msmrdMarkovModel msmrdMSM);

        void integrate(std::vector<particle> &parts) override;

        void computeTransitionsFromTransitionStates(std::vector<particle> &parts) override;

        void transition2BoundState(std::vector<particle> &parts, int iIndex, int jIndex, int endState) override;

        void applyEvents(std::vector<particle> &parts) override;


        /* Functions below might need to be overridden for more complex implementations. */

        //void computeTransitionsFromBoundStates(std::vector<particle> &parts) override;

        //void transition2UnboundState(std::vector<particle> &parts, int iIndex, int jIndex, int endState) override;

        //void transitionBetweenBoundStates(std::vector<particle> &parts, int iIndex, int jIndex, int endState) override;

        //void removeUnrealizedEvents(std::vector<particle> &parts) override;


        /* Functions exclusive to multiparticle MSM/RD*/

        bool acceptBindingEvent(std::vector<particle> &parts, int iIndex, int jIndex, int endState);

        void integrateDiffusionCompounds(std::vector<particle> &parts, double dt0);

        void addCompound(std::vector<particle> &parts, int iIndex, int jIndex, int endState);

        void updateParticlesInCompound(std::vector<particle> &parts, particleCompound &partCompound,
                vec3<double> deltar, quaternion<double> deltaq);

        void refreshParticlesInCompounds(std::vector<particle> &parts);

        void cleanParticleCompoundsVector(std::vector<particle> &parts);

        std::vector<int> findClosedBindingLoops(std::vector<particle> &parts);

        int getNumberOfBindingsInCompound(int compoundIndex);

        std::map<std::tuple<int,int>, int> getBindingsInCompound(int compoundIndex);

        int getCompoundSize(int compoundIndex);

    protected:

        /* Auxiliary functions used by functions above */

        void createCompound(std::vector<particle> &parts, int mainIndex, int secondIndex, int endState);

        void addParticleToCompound(std::vector<particle> &parts, int mainIndex, int secondIndex, int endState);

        void joinParticleCompounds(std::vector<particle> &parts, int mainIndex, int secondIndex, int endState);

        void alignParticlesInCompound(std::vector<particle> &parts, int mainIndex, int secondIndex);

        void closeCompound(std::vector<particle> &parts, int mainIndex, int secondIndex, int endState);


//        quaternion<double> getParticleCompoundOrientation(particle &mainParticle);

        /* Functions below might need to be overridden for more complex implementations. */

        //std::tuple<bool, bool> checkUnbindingCompounds(std::vector<particle> &parts, int iIndex, int jIndex);


    };

    template <typename templateMSM>
    msmrdMultiParticleIntegrator<templateMSM>::msmrdMultiParticleIntegrator(double dt, long seed,
            std::string particlesbodytype, int numParticleTypes,std::array<double,2> radialBounds,
            std::vector<templateMSM> MSMlist, msmrdMarkovModel msmrdMSM) :
            msmrdIntegrator<templateMSM>(dt, seed, particlesbodytype, numParticleTypes, radialBounds,
                  MSMlist, msmrdMSM) {
        discreteTrajClass = std::make_shared<patchyDimerTrajectory2> (patchyDimerTrajectory2(2,1));
            };

    template <typename templateMSM>
    msmrdMultiParticleIntegrator<templateMSM>::msmrdMultiParticleIntegrator(double dt, long seed,
            std::string particlesbodytype, int numParticleTypes, std::array<double,2> radialBounds,
            templateMSM MSMlist, msmrdMarkovModel msmrdMSM) :
            msmrdIntegrator<templateMSM>(dt, seed, particlesbodytype, numParticleTypes, radialBounds,
                                         MSMlist, msmrdMSM) {
        discreteTrajClass = std::make_shared<patchyDimerTrajectory2> (patchyDimerTrajectory2(2,1));
    };


    /**
    * Additional functions exclusive to multi-particle MSM/RD below
    */

    /* Integrates the diffusion of a compound of several particles. It assigns the new position and orientation
     * of the compound and also of the individual particles by calling the required functions (see template
     * functions in header). */
    template <typename templateMSM>
    void msmrdMultiParticleIntegrator<templateMSM>::integrateDiffusionCompounds(std::vector<particle> &parts, double dt0){
        for (auto &particleCompound : particleCompounds) {
            // Calculate change in position
            vec3<double> dr =  std::sqrt(2 * dt0 * particleCompound.D) * integrator::randg.normal3D(0, 1);
            // Calculate change in orientation
            vec3<double> dphi = std::sqrt(2 * dt0 * particleCompound.Drot) *
                    integrator::randg.normal3D(0, 1);
            quaternion<double> dquat = msmrdtools::axisangle2quaternion(dphi);
            // Update position and orientation of particles in compound
            updateParticlesInCompound(parts, particleCompound, dr,  dquat);
            // Update position and orientation of compound
            particleCompound.position += dr;
            particleCompound.orientation = dquat * particleCompound.orientation;
        }
    }


    /* Add particle complex into vector particleComplexes if particles bounded. If complex doesn't exist,
     * it creates it. If both particles belong to different complexes, it merges them. It also updates
     * the compoundIndex of each particle, see particle.hpp. It also assign the position and orientation
     * of the complex and the relative positions and orientations of all the particles within the complex. */
    template <typename templateMSM>
    void msmrdMultiParticleIntegrator<templateMSM>::addCompound(std::vector<particle> &parts, int iIndex,
                                                         int jIndex, int endState) {
        /* Assign indexes of particles depending if they belong to the main compound or the secondary compound. The
         * convention is the main compound is that with the smallest particle index between the binding partners. If no
         * compound exists, the mainIndex corresponds to the particle with the smallest index. If only one compound
         * exists the mainIndex is the one corresponding to the particle in the compound. */
        int mainIndex = 1 * iIndex;
        int secondIndex = 1 * jIndex;
        int endStateIndex = endState - 1;
        // If neither particle belongs to a compound, create one.
        if (parts[iIndex].compoundIndex == -1 and parts[jIndex].compoundIndex == -1) {
            createCompound(parts, mainIndex, secondIndex, endState);
        }
        // If one of the two doesn't belong to a compound, join the solo particle into the compound.
        else if (parts[iIndex].compoundIndex == -1 or parts[jIndex].compoundIndex == -1) {
            // Switch main and second index if neccesary (main particle must be in compound)
            if (parts[iIndex].compoundIndex == -1){
                // Inverts main and secondIndex and thus the def. of bound state (+1 to pass from index to boundstate)
                auto flippedBoundState = 1 + discreteTrajClass->getFlippedBoundStateIndex(endStateIndex);
                addParticleToCompound(parts, jIndex, iIndex, flippedBoundState);
            } else {
                addParticleToCompound(parts, iIndex, jIndex, endState);
            }
        }
        // If both belong to a compound, join the compounds together (keep the one with lower index)
        //else {
        else if (parts[iIndex].compoundIndex != parts[jIndex].compoundIndex){
                joinParticleCompounds(parts, mainIndex, secondIndex, endState);
            }
//        }
//        // If both belong to the same compound, close a loop in that complex
        else {
            closeCompound(parts, mainIndex, secondIndex, endState);
        }
         /* like the pentamer, where the binding of 5 particles together, essentially implies the pentamer formation. */
        // Deactivate corresponding particles since now their diffusion will be controlled by a particle compound */
        parts[mainIndex].deactivate();
        parts[secondIndex].deactivate();
    };

    /* Updates particles positions and orientations in a compound that diffused by deltar and rotated by deltaq */
    template <typename templateMSM>
    void msmrdMultiParticleIntegrator<templateMSM>::updateParticlesInCompound(std::vector<particle> &parts,
        particleCompound &partCompound, vec3<double> deltar, quaternion<double> deltaq) {
        // Update orientations of the particles. NOTE refParticle and compound always have same orientation!
        auto refOrientation = deltaq * partCompound.orientation;
        for (auto element : partCompound.relativeOrientations) {
            auto partIndex = element.first;
            auto relOrientation = element.second;
            parts[partIndex].orientation = refOrientation * relOrientation; //deltaOrientation;
        }
        // Update positions of the particles
        for (auto element : partCompound.relativePositions) {
            auto partIndex = element.first;
            auto relPosition = element.second;
            parts[partIndex].position = partCompound.position + deltar +
                                        msmrdtools::rotateVec(relPosition, deltaq * partCompound.orientation);
        }
    }

    /* Refreshes particles in all compounds. This means the nextPosition and nextOrientation are refreshed to avoid
     * conflicts wth th parent class. */
    template <typename templateMSM>
    void msmrdMultiParticleIntegrator<templateMSM>::refreshParticlesInCompounds(std::vector<particle> &parts) {
        // Update orientations of the particles. NOTE refParticle and compound always have same orientation!
        for (auto &partCompound : particleCompounds) {
            if (partCompound.isActive()) {
                for (auto element : partCompound.relativeOrientations) {
                    auto partIndex = element.first;
                    // Updates next orientation in case it is needed for parent functions calculations
                    parts[partIndex].setNextOrientation(parts[partIndex].orientation);
                }
                // Update positions of the particles
                for (auto element : partCompound.relativePositions) {
                    auto partIndex = element.first;
                    // Updates next position in case it is needed for parent functions calculations
                    parts[partIndex].setNextPosition(parts[partIndex].position);
                }
            }
        }
    }

    /* Deletes inactive complexes in particle complex vector, and updates indexes in particle list. Doesn't
     * need to do at every time step, but every now and then to free up memory. Not necessary but avoids possible
     * memory overflow in large systems. */
    template <typename templateMSM>
    void msmrdMultiParticleIntegrator<templateMSM>::cleanParticleCompoundsVector(std::vector<particle> &parts) {
        // Reverse iterator loop on particle compounds
        auto it = particleCompounds.end();
        while (it > particleCompounds.begin()) {
            it--;
            auto index = std::distance(particleCompounds.begin(), it);
            if (particleCompounds[index].active == false) {
                // Erase particle complex
                it = particleCompounds.erase(it);
                // Readjust the indexes of particles
                for (int j = 0; j < parts.size(); j++) {
                    // this should never happen
                    if (parts[j].compoundIndex == index) {
                        parts[j].compoundIndex = -1;
                    }
                    // move all the indices larger than i by -1 since we will delete one element
                    if (parts[j].compoundIndex > index) {
                        parts[j].compoundIndex--;
                    }
                }
            }
        }
    };

    /* Checks if there is any closed binding loop in any of the particle compounds. If so, it returns the size of
     * the loops found in a vector of intergers, which size is the number of loops */
    template <typename templateMSM>
    std::vector<int> msmrdMultiParticleIntegrator<templateMSM>::findClosedBindingLoops(std::vector<particle> &parts){
        std::vector<int> boundLoops = {};
        for (auto &particleCompound : particleCompounds) {
            if (particleCompound.active) {
                auto compoundSize = particleCompound.getSizeOfCompound();
                auto numBindings = particleCompound.getNumberOfbindings();
                if (compoundSize == numBindings) {
                    boundLoops.push_back(numBindings);
                }
                // Check for two-particle closed loop bindings (will not appear with method below )
                if (compoundSize == 2){
                    auto particleIndexes = particleCompound.boundPairsDictionary.begin()->first;
                    auto iIndex = std::get<0>(particleIndexes);
                    if (parts[iIndex].boundList[0] == parts[iIndex].boundList[1]) {
                        boundLoops.push_back(2);
                    }
                }
            }
        }
        return boundLoops;
    };

    /* Extracts number of bindings in a compound indexed by compoundIndex */
    template <typename templateMSM>
    int msmrdMultiParticleIntegrator<templateMSM>::getNumberOfBindingsInCompound(int compoundIndex) {
        return particleCompounds[compoundIndex].boundPairsDictionary.size();
    };

    /* Extracts bindings in a compound indexed by compoundIndex */
    template <typename templateMSM>
    std::map<std::tuple<int,int>, int> msmrdMultiParticleIntegrator<templateMSM>::getBindingsInCompound(int compoundIndex) {
        return particleCompounds[compoundIndex].boundPairsDictionary;
    };

    template <typename templateMSM>
    int msmrdMultiParticleIntegrator<templateMSM>::getCompoundSize(int compoundIndex) {
        return particleCompounds[compoundIndex].getSizeOfCompound();
    };




    /*
     * Auxiliary functions used by functions above
     */

    /* Creates a new compound from a binding between two particles. */
    template <typename templateMSM>
    void msmrdMultiParticleIntegrator<templateMSM>::createCompound(std::vector<particle> &parts, int mainIndex,
            int secondIndex, int endState) {
        int endStateIndex = endState - 1;
        particle &mainPart = parts[mainIndex];
        particle &secondPart = parts[secondIndex];
        std::tuple<int,int> pairIndices = std::make_tuple(mainIndex, secondIndex);
        std::map<std::tuple<int,int>, int> boundPairsDictionary = {{pairIndices, endState}};
        particleCompound pComplex = particleCompound(boundPairsDictionary);
        pComplex.referenceParticleIndex = mainIndex;
        /* Set relative positions and orientations in particle compound with respect to pentamer center, assuming
         * main particle is in origin with identity orientation. */
        pComplex.relativePositions.insert(std::pair<int, vec3<double>>(mainIndex, -1.0 * pentamerCenter));
        pComplex.relativeOrientations.insert(std::pair<int, quaternion<double>>(mainIndex,
                quaternion<double>(1,0,0,0)));
        // Get relative position and orientation with respect to main particle
        auto relPosition = discreteTrajClass->getRelativePosition(endStateIndex);
        auto relOrientation = discreteTrajClass->getRelativeOrientation(endStateIndex);
        // Set fixed position and orientation of newly bound particle
        secondPart.position = mainPart.position + msmrdtools::rotateVec(relPosition, mainPart.orientation);
        secondPart.orientation = mainPart.orientation * relOrientation;
        // Set particle compound position and orientation (as the center of a pentamer, same for all compounds)
        pComplex.position = mainPart.position + msmrdtools::rotateVec(pentamerCenter, mainPart.orientation);
        pComplex.orientation = mainPart.orientation; // make a drawing to understand why
        // Set relative position and orientation with respect to pentamer center (orientation w/respect to refParticle)
        relPosition += -1.0 * pentamerCenter;
        pComplex.relativePositions.insert(std::pair<int, vec3<double>>(secondIndex, relPosition));
        pComplex.relativeOrientations.insert(std::pair<int, quaternion<double>>(secondIndex, relOrientation));
        /* Set particle complex diffusion coefficients (This part could be extracted from MD data, here we simply
         * assume the same diffusion coefficient as the solo particles) */
        pComplex.D = 1.0 * mainPart.D;
        pComplex.Drot = 1.0 * mainPart.Drot;
        // Add particle compound to particle compound vector
        particleCompounds.push_back(pComplex);
        // Set new particle complex indices.
        mainPart.compoundIndex = static_cast<int>(particleCompounds.size() - 1);
        secondPart.compoundIndex = static_cast<int>(particleCompounds.size() - 1);
    };

    /* Adds a particle into an existing compound after a binding between a compound and a particle. The endState is
     * the bound state from the main particle (mainIndex) to the secondary particle (secondIndex)*/
    template <typename templateMSM>
    void msmrdMultiParticleIntegrator<templateMSM>::addParticleToCompound(std::vector<particle> &parts, int mainIndex,
            int secondIndex, int endState) {
        int endStateIndex = endState - 1;
        particle &mainPart = parts[mainIndex];
        particle &secondPart = parts[secondIndex];
        // Generate new binding description and insert it into compound.
        int compoundIndex = 1 * mainPart.compoundIndex;
        std::tuple<int,int> pairIndices = std::make_tuple(mainIndex, secondIndex);
        particleCompounds[compoundIndex].boundPairsDictionary.insert (
                std::pair<std::tuple<int,int>, int>(pairIndices, endState) );
        // Get relative position and orientation of new particle in complex
        auto relPosition = discreteTrajClass->getRelativePosition(endStateIndex);
        auto relOrientation = discreteTrajClass->getRelativeOrientation(endStateIndex);
        // Set fixed position and orientation of newly bound particle
        secondPart.position = mainPart.position + msmrdtools::rotateVec(relPosition, mainPart.orientation);
        secondPart.orientation = mainPart.orientation * relOrientation;
        // Set relative position w/respect to compound center and orientation w/respect to reference particle
        relPosition = particleCompounds[compoundIndex].relativePositions[mainIndex]
                + msmrdtools::rotateVec(relPosition, particleCompounds[compoundIndex].relativeOrientations[mainIndex]);
        relOrientation = particleCompounds[compoundIndex].relativeOrientations[mainIndex] * relOrientation;
        particleCompounds[compoundIndex].relativePositions.insert(std::pair<int, vec3<double>>(secondIndex,
                relPosition));
        particleCompounds[compoundIndex].relativeOrientations.insert(std::pair<int, quaternion<double>>(
                secondIndex, relOrientation));
        // Set new particle compound indices.
        secondPart.compoundIndex = 1 * compoundIndex;
    };

    /* Joins two compounds after a binding event. */
    template <typename templateMSM>
    void msmrdMultiParticleIntegrator<templateMSM>::joinParticleCompounds(std::vector<particle> &parts, int mainIndex,
            int secondIndex, int endState) {
        int endStateIndex = endState - 1;
        particle &mainPart = parts[mainIndex];
        particle &secondPart = parts[secondIndex];
        int mainCompoundIndex = 1 * mainPart.compoundIndex;
        int secondCompoundIndex = 1 * secondPart.compoundIndex;
        // Insert new binding into boundsPairs dicitionary
        std::tuple<int,int> pairIndices = std::make_tuple(mainIndex, secondIndex);
        particleCompounds[mainCompoundIndex].boundPairsDictionary.insert (
                std::pair<std::tuple<int,int>, int>(pairIndices, endState) );
        // Join complexes
        particleCompounds[mainCompoundIndex].joinCompound(particleCompounds[secondCompoundIndex]);
        // Set relative positions/orientations in particle compound, first for the particle where there was a binding.
        auto relPosition = discreteTrajClass->getRelativePosition(endStateIndex);
        auto relOrientation = discreteTrajClass->getRelativeOrientation(endStateIndex);
        // Set fixed position and orientation of newly bound particle
        secondPart.position = mainPart.position + msmrdtools::rotateVec(relPosition, mainPart.orientation);
        secondPart.orientation = mainPart.orientation * relOrientation;
        // Set relative position w/respect to compound center and orientation w/respect reference particle
        relPosition = particleCompounds[mainCompoundIndex].relativePositions[mainIndex]
                      + msmrdtools::rotateVec(relPosition,
                              particleCompounds[mainCompoundIndex].relativeOrientations[mainIndex]);
        relOrientation = particleCompounds[mainCompoundIndex].relativeOrientations[mainIndex] * relOrientation ;
        particleCompounds[mainCompoundIndex].relativePositions.insert(std::pair<int, vec3<double>>(secondIndex,
                                                                                               relPosition));
        particleCompounds[mainCompoundIndex].relativeOrientations.insert(std::pair<int, quaternion<double>>(
                secondIndex, relOrientation));
        /* ... then assign them to all the other particles in newly bound complex (aligns particles positions and
         * orientations of binding complex) */
        alignParticlesInCompound(parts, mainIndex, secondIndex);
        /* Deactivates secondary compound (empties dictionaries and flags complex to be deleted). */
        particleCompounds[secondCompoundIndex].deactivateCompound();
        // Set new particle complex indices to the one of the main compound index.
        secondPart.compoundIndex = 1 * mainCompoundIndex;
    };

    /* Aligns remaining particles in compound that just binded with another compound. This means rewriting positions
     * and orientations as well as relativePositions and Orientations in particleCompound. */
    template <typename templateMSM>
    void msmrdMultiParticleIntegrator<templateMSM>::alignParticlesInCompound(std::vector<particle> &parts,
                                                                             int mainIndex, int secondIndex) {
        particle &mainPart = parts[mainIndex];
        particle &secondPart = parts[secondIndex];
        int mainCompoundIndex = mainPart.compoundIndex;
        int secondCompoundIndex = secondPart.compoundIndex;
        std::map<std::tuple<int,int>, int> boundPairsDictionaryCopy(
                particleCompounds[secondCompoundIndex].boundPairsDictionary);
        int refParticleIndex = secondIndex;
        while (boundPairsDictionaryCopy.size() > 0) {
            for (auto it = boundPairsDictionaryCopy.cbegin(); it != boundPairsDictionaryCopy.cend(); it++) {
                auto index1 = std::get<0>(it->first);
                auto index2 = std::get<1>(it->first);
                auto state = it->second;
                auto stateIndex = state - 1;
                auto relPos = discreteTrajClass->getRelativePosition(stateIndex);
                auto relOrient = discreteTrajClass->getRelativeOrientation(stateIndex);
                if (index1 == refParticleIndex) {
                    parts[index2].position = parts[index1].position +
                            msmrdtools::rotateVec(relPos, parts[index1].orientation);
                    parts[index2].orientation = parts[index1].orientation * relOrient;
                    // Set relative positon/orientation
                    relPos = particleCompounds[mainCompoundIndex].relativePositions[index1]
                                  + msmrdtools::rotateVec(relPos,
                                          particleCompounds[mainCompoundIndex].relativeOrientations[index1]);
                    relOrient = particleCompounds[mainCompoundIndex].relativeOrientations[index1] * relOrient;
                    particleCompounds[mainCompoundIndex].relativePositions.insert(std::pair<int, vec3<double>>(
                            index2, relPos));
                    particleCompounds[mainCompoundIndex].relativeOrientations.insert(std::pair<int,
                            quaternion<double>>(index2, relOrient));
                    parts[index2].compoundIndex = mainCompoundIndex;
                    refParticleIndex = index2;
                    // Remove element from reference map and break for loop
                    it = boundPairsDictionaryCopy.erase(it);
                    break;
                }
                else if (index2 == refParticleIndex){
                    parts[index1].position = parts[index2].position +
                                             msmrdtools::rotateVec(relPos, parts[index2].orientation);
                    parts[index1].orientation = parts[index2].orientation * relOrient;
                    // Set relative positon/orientation
                    relPos = particleCompounds[mainCompoundIndex].relativePositions[index2]
                             + msmrdtools::rotateVec(relPos,
                                                     particleCompounds[mainCompoundIndex].relativeOrientations[index2]);
                    relOrient = particleCompounds[mainCompoundIndex].relativeOrientations[index2] * relOrient;
                    particleCompounds[mainCompoundIndex].relativePositions.insert(std::pair<int, vec3<double>>(
                            index1, relPos));
                    particleCompounds[mainCompoundIndex].relativeOrientations.insert(std::pair<int,
                            quaternion<double>>(index1, relOrient));
                    parts[index1].compoundIndex = mainCompoundIndex;
                    refParticleIndex = index1;
                    // Remove element from reference map and break loop
                    it = boundPairsDictionaryCopy.erase(it);
                    break;
                }
            }
        }
    };


    /* Closes a compound through a binding event that takes place within the same compound. Note relative positions
     * and orientations are not modified since this will in general be the end of the simulation. */
    template <typename templateMSM>
    void msmrdMultiParticleIntegrator<templateMSM>::closeCompound(std::vector<particle> &parts, int mainIndex,
                                                                          int secondIndex, int endState) {
        particle &mainPart = parts[mainIndex];
        particle &secondPart = parts[secondIndex];
        int endStateIndex = endState - 1;
        int mainCompoundIndex = mainPart.compoundIndex;
        // Insert new binding into boundsPairs dicitionary
        std::tuple<int,int> pairIndices = std::make_tuple(mainIndex, secondIndex);
        particleCompounds[mainCompoundIndex].boundPairsDictionary.insert (
                std::pair<std::tuple<int,int>, int>(pairIndices, endState) );
        // Set relative position and orientation
        // Set relative positions/orientations in particle compound, first for the particle where there was a binding.
        auto relPosition = discreteTrajClass->getRelativePosition(endStateIndex);
        auto relOrientation = discreteTrajClass->getRelativeOrientation(endStateIndex);
        // Set fixed position and orientation of newly bound particle
        secondPart.position = mainPart.position + msmrdtools::rotateVec(relPosition, mainPart.orientation);
        secondPart.orientation = mainPart.orientation * relOrientation;
        // Set relative position w/respect to compound center and orientation w/respect reference particle
        relPosition = particleCompounds[mainCompoundIndex].relativePositions[mainIndex]
                      + msmrdtools::rotateVec(relPosition,
                                              particleCompounds[mainCompoundIndex].relativeOrientations[mainIndex]);
        relOrientation = particleCompounds[mainCompoundIndex].relativeOrientations[mainIndex] * relOrientation ;
        particleCompounds[mainCompoundIndex].relativePositions.insert(std::pair<int, vec3<double>>(secondIndex,
                                                                                                   relPosition));
        particleCompounds[mainCompoundIndex].relativeOrientations.insert(std::pair<int, quaternion<double>>(
                secondIndex, relOrientation));
    };




//    /* Recovers the rotation quaternion for the particle compound based on the reference particle of the compound. It
//     * assumer the center of the compound is the center of a pentamer formation, so this is estremely dependent on the
//     * specific application. */
//    template <typename templateMSM>
//    quaternion<double> msmrdMultiParticleIntegrator<templateMSM>::getParticleCompoundOrientation(particle &mainParticle){
//        auto vec0 = vec3<double>(1.0/(2.0 * std::sin(M_PI/5.0)),0,0);
//        auto vec1 = vec3<double>(0,0,0);
//        auto vec2 = vec3<double>(0,1,0);
//        auto rotatedVec0 = msmrdtools::rotateVec(vec0,mainParticle.orientation);
//        auto rotatedVec1 = vec3<double>(0,0,0);
//        auto rotatedVec2 = msmrdtools::rotateVec(vec2,mainParticle.orientation);
//        auto outputQuat = msmrdtools::recoverRotationFromVectors(vec0,vec1,vec2,rotatedVec0,rotatedVec1,rotatedVec1);
//        return outputQuat;
//    }


//    /* Checks ...*/
//    std::tuple<bool, bool> msmrdMultiParticleIntegrator::checkUnbindingCompounds(std::vector<particle> &parts,
//                                                                               int iIndex, int jIndex) {
//        int cIndex = parts[iIndex].compoundIndex;
//
//        for (const auto &entry : particleCompounds[cIndex].boundPairsDictionary) {
//            auto key = entry.first;
//            auto value = entry.second
//        }
//
//
//    };


};

