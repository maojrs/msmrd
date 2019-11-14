//
// Created by maojrs on 11/5/19.
//

#include "integrators/msmrdPatchyProtein.hpp"

namespace msmrd {
    /**
     * Implementation of MSM/RD integrator for patchyProteins. Same constructor than msmrdIntegrator; however, it
     * calls a new discretization.
     */
    msmrdPatchyProtein::msmrdPatchyProtein(double dt, long seed, std::string particlesbodytype, int numParticleTypes,
                       std::array<double,2> radialBound, std::vector<ctmsm> MSMlist, msmrdMSM markovModel) :
            msmrdIntegrator<ctmsm>::msmrdIntegrator(dt, seed, particlesbodytype, numParticleTypes,
            radialBound, MSMlist, markovModel) {
        setPatchyProteinDiscretization();
    }

    msmrdPatchyProtein::msmrdPatchyProtein(double dt, long seed, std::string particlesbodytype, int numParticleTypes,
                       std::array<double,2> radialBound, ctmsm MSMlist, msmrdMSM markovModel) :
            msmrdIntegrator<ctmsm>::msmrdIntegrator(dt, seed, particlesbodytype, numParticleTypes,
                                                    radialBound, MSMlist, markovModel){
        setPatchyProteinDiscretization();
    };

    /* Calculates and set a new spherePartition and the positionOrientationParition of the msmrdIntegrator with
     * local parameters corresponding to the patchyProtein implementation default values. */
    void msmrdPatchyProtein::setPatchyProteinDiscretization() {
        int numSphericalSectionsPos = 6;
        int numRadialSectionsQuat = 4;
        int numSphericalSectionsQuat = 6;
        double relativeDistanceCutOff = radialBounds[1];

        auto tempPositionPart = std::make_shared<spherePartition>(spherePartition(numSphericalSectionsPos));
        auto tempPositionOrientationPart = std::make_shared<fullPartition> (fullPartition(relativeDistanceCutOff,
                numSphericalSectionsPos, numRadialSectionsQuat, numSphericalSectionsQuat));
        setDiscretization(tempPositionPart);
        setDiscretization(tempPositionOrientationPart);
    };

    /**
     * Auxiliary functions to be overriden
     */
    int msmrdPatchyProtein::setNewUnboundState(int unboundState, std::vector<particle> &parts, int partIndex) {
        int newState = 0;
        auto partType = parts[partIndex].type;
        if (MSMlist[partType].tmatrix.size() > 1) {
            parts[partIndex].activeMSM = true;
            newState = unboundState;
        }
        /* Set new unbound states (which also eliminates pair connections by resetting boundTo and
         * boundState to -1) and activate particles. */
        parts[partIndex].setState(newState);
        parts[partIndex].activate();
        // Sets diffusion coefficients of new unbound state
        auto diff = MSMlist[partType].Dlist[newState];
        auto diffRot = MSMlist[partType].Drotlist[newState];
        parts[partIndex].setDs(diff, diffRot);
        return newState;
    }


    /**
     * Main MSM/RD patchy protein integrator overridden functions
     */


    // Same as msmrdIntegrator<ctmsm>::computeCurrentTransitionState but adjusts output if part2.state == 1
    int msmrdPatchyProtein::computeCurrentTransitionState(particle &part1, particle &part2) {
        int currentTransitionState = msmrdIntegrator<ctmsm>::computeCurrentTransitionState(part1, part2);
        if (part2.state == 1 and currentTransitionState > markovModel.getMaxNumberBoundStates()) {
            if (rotation) {
                currentTransitionState += positionOrientationPart->numTotalSections;
            } else {
                currentTransitionState += positionPart->numSections;
            }
        }
        return currentTransitionState;
    }


    /* Same as msmrdIntegrator<ctmsm>::transition2UnboundState, with minor adjustments to deal with additional
     * transition states corresponding to part2.state=1. */
    void msmrdPatchyProtein::transition2UnboundState(std::vector<particle> &parts, int iIndex,
                                                     int jIndex, int endStateAlt) {
        int unboundState = 0;

        // Redefine endstate indexing, so it is understood by the partition/discretization.
        int index0 = markovModel.getMaxNumberBoundStates();
        int endState = endStateAlt - index0;

        /* If endState larger than sections in partition, take the module. Since transition states
         * with higher value correspond to same place in partition, but different state for particle2. Also
         * uses this to calculate unboundState of particle2. */
        if (endState > positionOrientationPart->numTotalSections) {
            unboundState = static_cast<int>(std::floor((endState - 1)/positionOrientationPart->numTotalSections));
            endState = (endState - 1) % positionOrientationPart->numTotalSections + 1;
        }

        // Calculates and sets next unbound states (of the unbound MSM). If no MSM, defaults to zero.
        setRandomUnboundState(parts, iIndex);
        setNewUnboundState(unboundState, parts, jIndex);

        // Extract relative position and orientation from partition and endstate
        auto relativePositionOrientation = getRelativePositionOrientation(endState);
        auto relPosition = std::get<0>(relativePositionOrientation);

        // Set next orientations based on the relative ones (parts[iIndex] keeps track of bound particle orientation)
        if (rotation) {
            auto relOrientation = std::get<1>(relativePositionOrientation);
            parts[iIndex].nextOrientation = 1.0 * parts[iIndex].orientation;
            parts[jIndex].nextOrientation = relOrientation * parts[iIndex].nextOrientation;
        }

        // Set next positions based on the relative ones (remember parts[iIndex] keeps track of bound particle position)
        parts[iIndex].nextPosition = parts[iIndex].position - 0.5*relPosition;
        parts[jIndex].nextPosition = parts[iIndex].nextPosition + relPosition;
    }

}