//
// Created by maojrs on 4/8/19.
//

#include "trajectories/discrete/patchyProteinTrajectory.hpp"

namespace msmrd {

    /*
     * Constructors and definition of bound states of patchy protein classes.
     * IMPORTANT: Need to define pointer to positionOrientationPartition and
     * call setBoundStates in constructors of childs of discreteTrajectory.
     */

    patchyProteinTrajectory::patchyProteinTrajectory(unsigned long Nparticles, int bufferSize) :
            discreteTrajectory(Nparticles, bufferSize) {

        setRadialBounds(1.25, 2.25);
        setTolerances(0.12, 0.12*2*M_PI);

        int numSphericalSectionsPos = 6; //7; //7;
        int numRadialSectionsQuat = 4; // 3; //5;
        int numSphericalSectionsQuat = 6; // 6; //7;
        positionOrientationPart = std::make_unique<positionOrientationPartition>(rUpperBound, numSphericalSectionsPos,
                                                                                 numRadialSectionsQuat,
                                                                                 numSphericalSectionsQuat);
        // Offset thetas in discretization to match patches
        double offTheta = M_PI/4.0;
        positionOrientationPart->setThetasOffset(offTheta);

        // Set bound states defined for the patchy protein example
        setBoundStates();
    };

    patchyProteinTrajectory::patchyProteinTrajectory(unsigned long Nparticles, int bufferSize,
                                                     double rLowerBound, double rUpperBound) :
            discreteTrajectory(Nparticles, bufferSize, rLowerBound, rUpperBound) {

        setTolerances(0.12, 0.12*2*M_PI);

        int numSphericalSectionsPos = 6; //7; //7;
        int numRadialSectionsQuat = 4; // 3; //5;
        int numSphericalSectionsQuat = 6; // 6; //7;
        positionOrientationPart = std::make_unique<positionOrientationPartition>(rUpperBound, numSphericalSectionsPos,
                                                                                 numRadialSectionsQuat,
                                                                                 numSphericalSectionsQuat);
        // Offset thetas in discretization to match patches
        double offTheta = M_PI/4.0;
        positionOrientationPart->setThetasOffset(offTheta);

        // Set bound states defined for the patchy protein example
        setBoundStates();
    };


    /* Sets bound states (metastable regions) of this patchy dimer implementation (two equal patcy particles,
    * each with two patches an angle angleDiff away with two stable relative orientations). The centers of the
    * metastable regions are given by a tuple of relative position and relative orientation. The size of the
    * regions are determined by tolerancePosition and toleranceOrientation*/
    void patchyProteinTrajectory::setBoundStates() {
        double angleDiff = 3 * M_PI / 5; // angle difference to form a pentamer
        /* Define relative position vectors from particle 1 at the origin. These two patches
         * point in the same direction as the two patches in the dimer. */
        std::array<vec3<double>, 6> relPos;
        relPos[0] = {1., 0., 0.};
        relPos[1] = {0., 1., 0.};
        relPos[2] = {0., 0., 1.};
        relPos[3] = {-1., 0., 0.};
        relPos[4] = {0., -1., 0.};
        relPos[5] = {0., 0., -1.};
        /* Relative rotations (assuming particle 1 fixed) of particle 2 that yield the 6 bound states
         * in the axis-angle representation. (One needs to make drawing to understand)*/
        std::array<vec3<double>, 6> rotations;
        rotations[0] = {0.0, 0.0, M_PI}; //ok
        rotations[1] = {0.0, 0.0, -M_PI / 2.0}; //ok
        rotations[2] = {0.0, M_PI / 2.0, 0.0}; //ok
        rotations[3] = {0.0, 0.0, 0.0}; //ok
        rotations[4] = {0.0, 0.0, M_PI / 2.0}; //ok
        rotations[5] = {0.0, -M_PI / 2.0, 0.0}; //ok
        /*Convert rotations in the axis angle representation to quaternions */
        std::array<quaternion<double>, 6> quatRotations;
        for (int i = 0; i < 6; i++) {
            quatRotations[i] = msmrdtools::axisangle2quaternion(rotations[i]);
        }
        /* Fill bound states with corresponding combinations of relative position vectors and quaternion orientations.
         * Note we need to take the conj, so it matches the relative orientation from particle 1, as used to define
         * the states in trajectories/discrete/discreteTrajectory.hpp */
        boundStates[0] = std::make_tuple(relPos[0], quatRotations[0].conj());
        boundStates[1] = std::make_tuple(relPos[1], quatRotations[1].conj());
        boundStates[2] = std::make_tuple(relPos[2], quatRotations[2].conj());
        boundStates[3] = std::make_tuple(relPos[3], quatRotations[3].conj());
        boundStates[4] = std::make_tuple(relPos[4], quatRotations[4].conj());
        boundStates[5] = std::make_tuple(relPos[5], quatRotations[5].conj());
    }


    /*
     * Alternative implementation of patchyProteinTrajectory
     */

    /* Based on parent class discreteTrajectory function. It differs with the original since it takes
     * the particle 2 state to choose a discrete state. It assumes particle can only bind, while particle 2
     * is in state 0. The previous implementation assumes the behavior of particle's 2 state is averaged by
     * the MSM. */
    int patchyProteinTrajectory2::sampleDiscreteState(particle part1, particle part2) {
        // Initialize sample with value zero (unbound state)
        int discreteState = 0;

        /* Calculate relative position taking into account periodic boundary measured
         * from i to j (gets you from i to j). */
        vec3<double> relativePosition = calculateRelativePosition(part1.position, part2.position);

        // Rotate relative position to match the reference orientation of particle 1. (VERY IMPORTANT)
        relativePosition = msmrdtools::rotateVec(relativePosition, part1.orientation.conj());
        quaternion<double> quatReference = {1,0,0,0}; // we can then define reference quaternion as identity.

        // Calculate relative orientation (w/respect to particle 1)
        quaternion<double> relativeOrientation;
        //relativeOrientation = part1.orientation.conj() * part2.orientation;
        relativeOrientation =  part2.orientation * part1.orientation.conj();


        // Extract current state, save into sample and return sample
        int secNum;
        if (relativePosition.norm() < rLowerBound) {
            // Only sample bound states if part2 is in state 0.
            if (part2.state == 0) {
                discreteState = getBoundState(relativePosition, relativeOrientation);
            } else{
                // Returns -1 so functions in discretizeTrajectory can usbstitute with prevsample if using coreMSM.
                discreteState = -1;
            }
        } else if (relativePosition.norm() < positionOrientationPart->relativeDistanceCutOff) {
            // Get corresponding section numbers from spherical partition to classify its state
            secNum = positionOrientationPart->getSectionNumber(relativePosition, relativeOrientation, quatReference);
            // Take into account the state of particle 2 to define state numbering
            secNum += part2.state * positionOrientationPart->numTotalSections;
            // Make sure bound states and transitions states correspond to different numbers
            discreteState  = maxNumberBoundStates + secNum;
        }
        return discreteState;
    }

}
