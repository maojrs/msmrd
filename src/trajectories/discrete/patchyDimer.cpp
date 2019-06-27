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
        spherePart = std::make_unique<spherePartition>(numSphericalSectionsPos);
        quaternionPart = std::make_unique<quaternionPartition>(numRadialSectionsQuat, numSphericalSectionsPos);
        positionOrientationPart = std::make_unique<positionOrientationPartition>(radialCutOff, numSphericalSectionsPos,
                                                                                 numRadialSectionsQuat,
                                                                                 numSphericalSectionsQuat);
        angularStates = numSphericalSectionsPos*(numSphericalSectionsPos + 1)/2;
        setMetastableRegions();
    };


    void patchyDimer::sampleDiscreteTrajectory(double time, std::vector<particle> &particleList) {
        // Initialize sample with value zero
        std::vector<int> sample{0};

        // Relative position and orientation measured from i to j (gets you from i to j)
        vec3<double> relativePosition;
        quaternion<double> relativeOrientation;

        // Section number in positionOrientationPartition
        int secNum;

        // Extract discrete trajectory from 2 particle simulation, i = particle 1, j = particle 2
        int i = 0;
        int j = 1;

        // Calculate relative orientation
        auto orientation1 = particleList[i].orientation;
        auto orientation2 = particleList[j].orientation;
        relativeOrientation = orientation2*orientation1.conj();

        // Calculate relative distance taking into account periodic boundary.
        relativePosition = msmrdtools::calculateRelativePosition(particleList[i].position, particleList[j].position,
                boundaryActive, domainBoundary->getBoundaryType(), domainBoundary->boxsize);

        // Extract current state, save into sample and push sample to discreteTrajectoryData.
        if (relativePosition.norm() < 1.25) {
            sample = std::vector<int>{ getBoundState(orientation1, orientation2) };
//            if (sample[0] == -1) {
//                secNum = positionOrientationPart->getSectionNumber(relativePosition, relativeOrientation, orientation1);
//                sample  = std::vector<int>{ maxNumberBoundStates + secNum };
//            }
        } else if (relativePosition.norm() < positionOrientationPart->relativeDistanceCutOff) {
            // Get corresponding section numbers from spherical partition to classify its state
            secNum = positionOrientationPart->getSectionNumber(relativePosition, relativeOrientation, orientation1);
            sample  = std::vector<int>{maxNumberBoundStates + secNum};
        }
        prevsample = sample[0];
        discreteTrajectoryData.push_back(sample);
    };


    /* Given two orientations, return either bound state A (1) or bound state B (2).
     * Note it can be modified to take into account symmetries, using the symmetryQuaternions. */
    int patchyDimer::getBoundState(quaternion<double> q1, quaternion<double> q2) {
        std::array<quaternion<double>,8> relOrientations;
        // Calculate equivalent relative orientations, measured from particle 1
        relOrientations[0] = q2 * q1.conj();
        relOrientations[1] = (q2*symmetryQuaternions[0]) * q1.conj();
        relOrientations[2] = q2 * (q1*symmetryQuaternions[0]).conj();
        relOrientations[3] = (q2*symmetryQuaternions[0]) * (q1*symmetryQuaternions[0]).conj();
        // Calculate equivalent relative orientations, measured from particle 2
        relOrientations[4] = relOrientations[0].conj();
        relOrientations[5] = relOrientations[1].conj();
        relOrientations[6] = relOrientations[2].conj();
        relOrientations[7] = relOrientations[3].conj();

        // Looping over all equivalent relative orientations, determines if it is in state A (1) or B (2)
        for (auto &relOrient : relOrientations) {
            if (msmrdtools::quaternionAngleDistance(relOrient, rotMetastableStates[0]) < 2*M_PI*tolerance) {
                return 1;
            }
            if (msmrdtools::quaternionAngleDistance(relOrient, rotMetastableStates[1]) < 2*M_PI*tolerance) {
                return 2;
            }
        }
        return prevsample; //-1;
    };


    /* Given two position and orientations, return either bound state A (1), bound state B (2) or unbound (-1).
     * This is a special version for PyBind, useful when calculating benchmarks in python interface. */
    int patchyDimer::getBoundStatePyBind(particle part1, particle part2) {

        // Calculate relative orientation
        vec3<double> relativePosition;
        auto position1 = part1.position;
        auto position2 = part2.position;
        auto orientation1quat = part1.orientation;
        auto orientation2quat = part2.orientation;

        // Calculate relative distance taking into account periodic boundary.
        if (boundaryActive) {
            relativePosition = msmrdtools::calculateRelativePosition(position1, position2,
                                                                                  boundaryActive,
                                                                                  domainBoundary->getBoundaryType(),
                                                                                  domainBoundary->boxsize);
        } else {
            relativePosition = position2 - position1;
        }

        // Extract current state, save into sample and push sample to discreteTrajectoryData.
        int boundState = -1;
        if (relativePosition.norm() < 1.25) {
            boundState = getBoundState(orientation1quat, orientation2quat);
        }
        if (boundState != 1 and boundState != 2) {
            boundState = -1;
        }
        return boundState;
    };

    /*
     * Define metastable relative orientations (including symmetric quaternion)
     */
    void patchyDimer::setMetastableRegions() {
        // Define symmetry quaternions that denotes rotational symmetry of particle, if needed
        //symmetryQuaternions.resize(1); //resize if needed
        symmetryQuaternions[0] = msmrdtools::axisangle2quaternion(vec3<double>(M_PI, 0, 0));
        // Define 2 main relative rotations (quaternions) that define bounds states A and B.
        const int numStates = 2;
        std::array<vec3<double>, numStates> axisAngleVecs;
        rotMetastableStates.resize(numStates);
        axisAngleVecs[0] = vec3<double>(0, 0, M_PI - 3.0*M_PI/5);
        axisAngleVecs[1] = vec3<double>(0, 0, M_PI);
        for (int i = 0; i < numStates; i++){
            rotMetastableStates[i] = msmrdtools::axisangle2quaternion(axisAngleVecs[i]);
        }
    }

}