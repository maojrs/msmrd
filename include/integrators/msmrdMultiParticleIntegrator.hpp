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
        /**
         * @param particleCompounds: vector containing particle compounds (two or more particles bound
         * together) to track them and diffuse them.
         * @param discreteTrajClass: pointer to trajectory class used for the discretization. We need this to extract
         * the bound states and their relative positions and orientations. The specific discrete trajectory
         * class will be fixed in the constructor, so this needs to be modified for different MSM/RD applications.
         */

        msmrdMultiParticleIntegrator(double dt, long seed, std::string particlesbodytype, int numParticleTypes,
                        std::array<double,2> radialBounds, std::vector<templateMSM> MSMlist,
                        msmrdMarkovModel msmrdMSM);

        msmrdMultiParticleIntegrator(double dt, long seed, std::string particlesbodytype, int numParticleTypes,
                        std::array<double,2> radialBounds, templateMSM MSMlist, msmrdMarkovModel msmrdMSM);

        using msmrdIntegrator<ctmsm>::msmrdIntegrator;

        void integrate(std::vector<particle> &parts) override;

        void computeTransitionsFromTransitionStates(std::vector<particle> &parts) override;

        void computeTransitionsFromBoundStates(std::vector<particle> &parts) override;

        void transition2BoundState(std::vector<particle> &parts, int iIndex, int jIndex, int endState) override;

        void transition2UnboundState(std::vector<particle> &parts, int iIndex, int jIndex, int endState) override;

        //void transitionBetweenBoundStates(std::vector<particle> &parts, int iIndex,
        //                                  int jIndex, int endState) override;

        //void removeUnrealizedEvents(std::vector<particle> &parts) override;

    protected:

        /* Functions exclusive to multiparticle MSM/RD*/

        void integrateDiffusionCompounds(double dt0, std::vector<particle> &parts) {};

        int addCompound(std::vector<particle> &parts, int iIndex, int jIndex, int endState);

        void setCompoundPositionOrientation(std::vector<particle> &parts, int iIndex, int jIndex,
                                            int mainComplexSize);

        void setParticlePositionOrientation(std::vector<particle> &parts, particleCompound &particleCompound,
                vec3<double> deltar,quaternion<double> deltaq, int particleIndex);





//        std::tuple<bool, bool> checkUnbindingCompounds(std::vector<particle> &parts, int iIndex, int jIndex);

        void updateParticleComplexesVector(std::vector<particle> &parts);

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


    /* Add particle complex into vector particleComplexes if particles bounded. If complex doesn't exist,
     * it creates it. If both particles belong to different complexes, it merges them. It also updates
     * the compoundIndex of each particle, see particle.hpp. It also returns an integer corresponding to
     * the size of the main compound in the case of two compounds joining, or simply one otherwise. */
    template <typename templateMSM>
    int msmrdMultiParticleIntegrator<templateMSM>::addCompound(std::vector<particle> &parts, int iIndex,
                                                         int jIndex, int endState) {
        int mainCompoundSize = 1;
        particle &iPart = parts[iIndex];
        particle &jPart = parts[jIndex];
        // If neither particle belongs to a complex, create one
        if (iPart.compoundIndex == -1 and jPart.compoundIndex == -1) {
            std::tuple<int,int> pairIndices = std::make_tuple(iIndex, jIndex);
            std::map<std::tuple<int,int>, int> boundPairsDictionary = {{pairIndices, endState}};
            particleCompound pComplex = particleCompound(boundPairsDictionary);
            particleCompounds.push_back(pComplex);
            // Set new particle complex indices.
            iPart.compoundIndex = static_cast<int>(particleCompounds.size() - 1);
            jPart.compoundIndex = static_cast<int>(particleCompounds.size() - 1);
        }
            // If one of the two doesn't belong to a complex, join the solo particle into the complex.
        else if (iPart.compoundIndex * jPart.compoundIndex < 0) {
            // Generate new binding description
            int compoundIndex = std::max(iPart.compoundIndex, jPart.compoundIndex);
            std::tuple<int,int> pairIndices = std::make_tuple(iIndex, jIndex);
            //Insert new binding decription into complex
            particleCompounds[compoundIndex].boundPairsDictionary.insert (
                    std::pair<std::tuple<int,int>, int>(pairIndices, endState) );
            // Extract main compound size and increase complex size by one.
            mainCompoundSize = particleCompounds[compoundIndex].compoundSize;
            particleCompounds[compoundIndex].compoundSize ++;
            // Set new particle complex indices.
            iPart.compoundIndex = compoundIndex;
            jPart.compoundIndex = compoundIndex;
        }
            // If both belong to a complex, join the complexes together (keep the one with lower index)
        else {
            // Add new bound pair to particle complex
            int iCompoundIndex = std::min(iPart.compoundIndex, jPart.compoundIndex);
            int jCompoundIndex = std::max(iPart.compoundIndex, jPart.compoundIndex);
            std::tuple<int,int> pairIndices = std::make_tuple(iIndex, jIndex);
            particleCompounds[iCompoundIndex].boundPairsDictionary.insert (
                    std::pair<std::tuple<int,int>, int>(pairIndices, endState) );
            // Join complexes and flag complex with larger index to be deleted.
            particleCompounds[iCompoundIndex].joinParticleCompound(particleCompounds[jCompoundIndex]);
            particleCompounds[jCompoundIndex].active = false;
            // Extract main compound size, increase complex size by one and make zero the empty one
            mainCompoundSize = particleCompounds[iCompoundIndex].compoundSize;
            particleCompounds[iCompoundIndex].compoundSize ++;
            particleCompounds[jCompoundIndex].compoundSize = 0;
            // Set new particle complex indices to the one of the lower index.
            iPart.compoundIndex = iCompoundIndex;
            jPart.compoundIndex = iCompoundIndex;
        }
        return mainCompoundSize;
    };


    /* Sets compound position to the average, its orientation is always set up initially at one and only changed
     * by diffusion. It also sets the position reference and orientation reference for a newly created compound. */
    template <typename templateMSM>
    void msmrdMultiParticleIntegrator<templateMSM>::setCompoundPositionOrientation(std::vector<particle> &parts,
            int iIndex, int jIndex, int mainCompoundSize) {
        int compoundIndex = parts[iIndex].compoundIndex;
        int compoundSize = particleCompounds[compoundIndex].compoundSize;
        if (compoundSize == 2) {
            // Set average position
            particleCompounds[compoundIndex].position = 0.5*(parts[iIndex].position + parts[jIndex].position);
            //vec3<double> diffVec = parts[jIndex].position - parts[iIndex].position;
            particleCompounds[compoundIndex].positionReference = parts[iIndex].position; // - 0.5*diffVec/diffVec.norm();
            particleCompounds[compoundIndex].orientationReference = 1.0*parts[iIndex].orientation;
        } else {
            particleCompounds[compoundIndex].position = (mainCompoundSize * particleCompounds[compoundIndex].position +
                                                         (compoundSize - mainCompoundSize) *
                                                         parts[jIndex].position)/compoundSize;
            particleCompounds[compoundIndex].positionReference = parts[iIndex].position; // - 0.5*diffVec/diffVec.norm();
            particleCompounds[compoundIndex].orientationReference = 1.0*parts[iIndex].orientation;
        }
    };

    template <typename templateMSM>
    void msmrdMultiParticleIntegrator<templateMSM>::setParticlePositionOrientation(std::vector<particle> &parts,
            particleCompound &particleCompound, vec3<double> deltar, quaternion<double> deltaq, int particleIndex) {
        if (particleIndex == particleCompound.referenceIndex) {
            auto vec0 = parts[particleIndex].position;
            auto offAxisPoint = particleCompound.position;
            auto rotatedVec0 = deltar + msmrdtools::rotateVecOffAxis(vec0, deltaq, offAxisPoint);
            parts[particleIndex].position = rotatedVec0;
            auto vec1 = vec3<double>(1, 1, 0); //reference vector 1 to control orientation
            auto vec2 = vec3<double>(1, -1, 0); //reference vector 2 to control orientation
            vec1 = msmrdtools::rotateVec(vec1, parts[particleIndex].orientation);
            vec2 = msmrdtools::rotateVec(vec2, parts[particleIndex].orientation);
            auto rotatedVec1 = deltar + msmrdtools::rotateVecOffAxis(vec1, deltaq, offAxisPoint);
            auto rotatedVec2 = deltar + msmrdtools::rotateVecOffAxis(vec2, deltaq, offAxisPoint);
            quaternion<double> newOrientation = msmrdtools::recoverRotationFromVectors(vec0, vec1, vec2,
                                                                                       rotatedVec0, rotatedVec1, rotatedVec2);
            parts[particleIndex].orientation = newOrientation;
        }
    }


    /* Deletes inactive complexes in particle complex vector, and updates indexes in particle list. Doesn't
     * need to do at every time step, but every now and then to free up memory. */
    template <typename templateMSM>
    void msmrdMultiParticleIntegrator<templateMSM>::updateParticleComplexesVector(std::vector<particle> &parts) {
        for (size_t i = particleCompounds.size(); i--;) {
            if (particleCompounds[i].active == false) {
                // Erase particle complex
                particleCompounds.erase(particleCompounds.begin() + i);
                // Readjust the indexes of particles
                for (int j = 0; j < parts.size(); j++) {
                    // this should never happen
                    if (parts[j].compoundIndex == i) {
                        parts[j].compoundIndex = -1;
                    }
                    // move all the indices larger than i by -1 since we will delete one element
                    if (parts[j].compoundIndex > i) {
                        parts[j].compoundIndex--;
                    }
                }
            }
        }
    };


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

