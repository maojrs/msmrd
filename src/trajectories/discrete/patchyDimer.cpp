//
// Created by maojrs on 1/22/19.
//
#include "trajectories/discrete/patchyDimer.hpp"

namespace msmrd {

    /*
     * Chooses how to discretize the full trajectory (trajectoryData) into a discretized trajectory
     * to be analyzed and extracted into a Markov state model. This specific discretization will follow
     * the core MSM approach
     */
    patchyDimer::patchyDimer(unsigned long Nparticles, int bufferSize) :
    trajectoryPositionOrientation(Nparticles, bufferSize) {
        double radialCutOff = 2.2;
        int numSphericalSectionsPos = 7;
        int numRadialSectionsQuat = 5;
        int numSphericalSectionsQuat = numSphericalSectionsPos;
        discreteTrajectoryData.reserve(bufferSize);
        positionOrientationPart = std::make_unique<positionOrientationPartition>(radialCutOff, numSphericalSectionsPos,
                                                                                 numRadialSectionsQuat,
                                                                                 numSphericalSectionsQuat);
        setMetastableRegions();
    };


    void patchyDimer::sampleDiscreteTrajectory(double time, std::vector<particle> &particleList) {
        // Initialize sample with value zero
        std::vector<int> sample{0};

        auto part1 = particleList[0];
        auto part2 = particleList[1];

        /* Calculate relative position taking into account periodic boundary measured
         * from i to j (gets you from i to j). */
        vec3<double> relativePosition;
        if (boundaryActive) {
            relativePosition = msmrdtools::calculateRelativePosition(part1.position, part2.position,
                    boundaryActive, domainBoundary->getBoundaryType(), domainBoundary->boxsize);
        } else {
            relativePosition = part2.position - part1.position;
        }
        // Rotate relative position to match the reference orientation of particle 1. (VERY IMPORTANT)
        relativePosition = msmrdtools::rotateVec(relativePosition, part1.orientation.conj());
        quaternion<double> quatReference = {1,0,0,0}; // we can then define reference quaternion as identity.

        // Calculate relative orientation
        quaternion<double> relativeOrientation;
        relativeOrientation = part1.orientation.conj() * part2.orientation;

        // Extract current state, save into sample and push sample to discreteTrajectoryData.
        int secNum;
        if (relativePosition.norm() < 1.25) {
            sample = std::vector<int>{ getBoundState(relativePosition, relativeOrientation) };
        } else if (relativePosition.norm() < positionOrientationPart->relativeDistanceCutOff) {
            // Get corresponding section numbers from spherical partition to classify its state
            secNum = positionOrientationPart->getSectionNumber(relativePosition, relativeOrientation, quatReference);
            sample  = std::vector<int>{maxNumberBoundStates + secNum};
        }
        prevsample = sample[0];
        discreteTrajectoryData.push_back(sample);
    };


    /* Given two particles, use their positions and orientations to determine if they are in one of
     * the 8 bound states (1 to 8). If not, return the value of the previous state */
    int patchyDimer::getBoundState(vec3<double> relativePosition, quaternion<double> relativeOrientation) {

        /* Check if it matches a bound states, if so return the corresponding state. Otherwise
         * return the previous state. */
        vec3<double> relPosCenter;
        quaternion<double> relQuatCenter;
        double angleDistance;
        for (int i = 0; i < 8; i++) {
            relPosCenter = std::get<0>(boundStates[i]);
            relQuatCenter = std::get<1>(boundStates[i]);

            if ( (relPosCenter - relativePosition).norm() <= tolerancePosition) {
                angleDistance = msmrdtools::quaternionAngleDistance(relQuatCenter, relativeOrientation);
                if  ( angleDistance < toleranceOrientation) {
                    return i + 1;
                }
            }
        }

        return prevsample;
    };


    /* Similar to getBoundState, but return -1 when they don't correspond to any bound state. This is
     * a special version for PyBind, useful when calculating benchmarks in python interface. Note it is a mixture
     * between sampleDiscreteTrjacteory and getBoundState. */
    int patchyDimer::getBoundStatePyBind(particle part1, particle part2) {
        // Calculate relative distance taking into account periodic boundary.
        vec3<double> relativePosition;
        if (boundaryActive) {
            relativePosition = msmrdtools::calculateRelativePosition(part1.position, part2.position,
                    boundaryActive, domainBoundary->getBoundaryType(), domainBoundary->boxsize);
        } else {
            relativePosition = part2.position - part1.position;
        }
        // Rotate relative position to match the reference orientation of particle 1.
        relativePosition = msmrdtools::rotateVec(relativePosition, part1.orientation.conj());

        // Calculate relative orientation
        quaternion<double> relativeOrientation = part1.orientation.conj() * part2.orientation;

        // Check if it matches a bound state, if so return the corresponding state. Otherwise return -1.
        vec3<double> relPosCenter;
        quaternion<double> relQuatCenter;
        double angleDistance;
        for (int i = 0; i < 8; i++) {
            relPosCenter = std::get<0>(boundStates[i]);
            relQuatCenter = std::get<1>(boundStates[i]);

            if ( (relPosCenter - relativePosition).norm() <= tolerancePosition) {
                angleDistance = msmrdtools::quaternionAngleDistance(relQuatCenter, relativeOrientation);
                if  ( angleDistance < toleranceOrientation) {
                    return i + 1;
                }
            }
        }

        return -1;
    };


    void patchyDimer::setMetastableRegions() {
        double angleDiff = 3 * M_PI / 5; // angle difference to form a pentamer
        /* Define relative position vectors from particle 1 at the origin. These two patches
         * point in the same direction as the two patches in the dimer. */
        vec3<double> relPos1 = {std::cos(angleDiff / 2.0), std::sin(angleDiff / 2.0), 0};
        vec3<double> relPos2 = {std::cos(angleDiff / 2.0), std::sin(-angleDiff / 2.0), 0};
        vec3<double> relPos1orthogonal = {-1.0 * std::sin(angleDiff / 2.0), std::cos(angleDiff / 2.0), 0.0};
        vec3<double> relPos2orthogonal = {std::sin(angleDiff / 2.0), std::cos(angleDiff / 2.0), 0.0};
        /* Relative rotations (from particle 1) of particle 2 that yield the 8 bound states
         * in the axis-angle representation. (One needs to make drawing to understand)*/
        std::array<vec3<double>, 8> rotations;
        rotations[0] = M_PI * relPos1orthogonal; //ok
        rotations[1] = {0.0, 0.0, -2 * M_PI / 5.0}; //ok
        rotations[2] = {0.0, 0.0, M_PI}; //ok
        rotations[3] = {0.0, M_PI, 0.0}; //ok
        // --first 4 rotations correspond to binding on top patch of particle 1, next 4 rotations to bottom patch
        rotations[4] = M_PI * relPos2orthogonal; //ok
        rotations[5] = {0.0, 0.0, 2 * M_PI / 5.0}; //ok
        rotations[6] = {0.0, 0.0, M_PI}; //ok
        rotations[7] = {0.0, M_PI, 0.0}; //ok
        /*Convert rotations in the axis angle representation to quaternions */
        std::array<quaternion<double>, 8> quatRotations;
        for (int i = 0; i < 8; i++) {
            quatRotations[i] = msmrdtools::axisangle2quaternion(rotations[i]);
        }
        // Fill bound states with corresponding combinations of relative position vectors and quaternion orientations.
        boundStates[0] = std::make_tuple(relPos1, quatRotations[0]);
        boundStates[1] = std::make_tuple(relPos1, quatRotations[1]);
        boundStates[2] = std::make_tuple(relPos1, quatRotations[2]);
        boundStates[3] = std::make_tuple(relPos1, quatRotations[3]);
        boundStates[4] = std::make_tuple(relPos2, quatRotations[4]);
        boundStates[5] = std::make_tuple(relPos2, quatRotations[5]);
        boundStates[6] = std::make_tuple(relPos2, quatRotations[6]);
        boundStates[7] = std::make_tuple(relPos2, quatRotations[7]);
    }

}