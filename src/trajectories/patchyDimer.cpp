//
// Created by maojrs on 1/22/19.
//
#include "trajectories/patchyDimer.hpp"

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

        // Calculate relative distance. If particle within periodic boundary, take that into account.
        if (boundaryActive and domainBoundary->getBoundaryType() == "periodic") {
            auto boxsize = domainBoundary->boxsize;
            relativePosition = msmrdtools::distancePeriodicBox(particleList[j].position, particleList[i].position, boxsize);
        } else {
            relativePosition = particleList[j].position - particleList[i].position;
        }

        // Extract current state, save into sample and push sample to discreteTrajectoryData.
        if (relativePosition.norm() < 1.25) {
            sample = std::vector<int>{ getBoundState(orientation1, orientation2) };
//            if (sample[0] == -1) {
//                secNum = positionOrientationPart->getSectionNumber(relativePosition, relativeOrientation, orientation1);
//                sample  = std::vector<int>{ 10 + secNum };
//            }
        } else if (relativePosition.norm() < positionOrientationPart->relativeDistanceCutOff) {
            // Get corresponding section numbers from spherical partition to classify its state
            secNum = positionOrientationPart->getSectionNumber(relativePosition, relativeOrientation, orientation1);
            sample  = std::vector<int>{10 + secNum};
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

//    /* Given two orientations, return either bound state A (1) or bound state B (2).
//     * Note it can be modified to take into account symmetries, using the symmetryQuaternions. */
//    int patchyDimer::getBoundState(quaternion<double> q1, quaternion<double> q2) {
//        std::array<quaternion<double>,8> relOrientations;
//        // Calculate equivalent relative orientations, measured from particle 1
//        relOrientations[0] = q2 * q1.conj();
//        relOrientations[1] = (q2*symmetryQuaternions[0]) * q1.conj();
//        relOrientations[2] = q2 * (q1*symmetryQuaternions[0]).conj();
//        relOrientations[3] = (q2*symmetryQuaternions[0]) * (q1*symmetryQuaternions[0]).conj();
//        // Calculate equivalent relative orientations, measured from particle 2
//        relOrientations[4] = relOrientations[0].conj();
//        relOrientations[5] = relOrientations[1].conj();
//        relOrientations[6] = relOrientations[2].conj();
//        relOrientations[7] = relOrientations[3].conj();
//
//        // Looping over all equivalent relative orientations, determines if it is in state A (1) or B (2)
//        for (auto &relOrient : relOrientations) {
//            if (msmrdtools::quaternionAngleDistance(relOrient, rotMetastableStates[0]) < 2*M_PI*tolerance) {
//                return 1;
//            }
//            if (msmrdtools::quaternionAngleDistance(relOrient, rotMetastableStates[1]) < 2*M_PI*tolerance) {
//                return 2;
//            }
//        }
//        return prevsample; //-1;
//    };

//    void patchyDimer::sampleDiscreteTrajectory(double time, std::vector<particle> &particleList) {
//        std::vector<int> sample;
//        vec3<double> relativePosition; // Measured from i to j
//        quaternion<double> orientation1;
//        quaternion<double> orientation2;
//        //quaternion<double> relativeOrientation; // Measured from i to j
//        // Rotated relative position vectors w/respect to each particle
//        vec3<double> rotatedRij;
//        vec3<double> rotatedRji;
//        // Section numbers corresponding to the relative orientations in the sphere Partition
//        int secNum1;
//        int secNum2;
//
//        // Extract discrete trajectory from 2 particle simulation
//        int i = 0; // index of particle 1
//        int j = 1; // index of particle 2
//        sample = std::vector<int>{0};
//        // Extract info from two rlevant particles in list
//        orientation1 = particleList[i].orientation;
//        orientation2 = particleList[j].orientation;
//        // Calculate relative distance. If box periodic boundary, take that into account.
//        if (boundaryActive and domainBoundary->getBoundaryType() == "periodic") {
//            auto boxsize = domainBoundary->boxsize;
//            relativePosition = msmrdtools::distancePeriodicBox(particleList[j].position, particleList[i].position, boxsize);
//        } else {
//            relativePosition = particleList[j].position - particleList[i].position;
//        }
//        // Rotate the relative position vector to the frame of reference fixed in each of the two particles.
//        rotatedRij = msmrdtools::rotateVec(relativePosition, orientation1);
//        rotatedRji = msmrdtools::rotateVec(-1*relativePosition, orientation2);
//        // Evaluate if the two particles are in a bound state
//        if (relativePosition.norm() < 1.2) {
//            sample = std::vector<int>{ getBoundState(orientation1, orientation2) };
//            if (sample[0] == -1) {
//                secNum1 = spherePart->getSectionNumber(rotatedRij);
//                secNum2 = spherePart->getSectionNumber(rotatedRji);
//                sample  = std::vector<int>{ 1000 + std::min(secNum1,secNum2)*10 + std::max(secNum1,secNum2) };
//            }
//        } else if (relativePosition.norm() < 2.0) {
//            // Get corresponding section numbers from spherical partition to classify its state
//            secNum1 = spherePart->getSectionNumber(rotatedRij);
//            secNum2 = spherePart->getSectionNumber(rotatedRji);
//            sample  = std::vector<int>{ std::min(secNum1,secNum2)*10 + std::max(secNum1,secNum2) };
//        }
//        prevsample = sample[0];
//        discreteTrajectoryData.push_back(sample);
//
//    };
//
//    /* Given two orientations, return either bound state A (1) or bound state B (2).
// * Note it can be modified to take into account symmetries, using the symmetryQuaternions. */
//    int patchyDimer::getBoundState(quaternion<double> q1, quaternion<double> q2) {
//        std::array<quaternion<double>,8> relOrientations;
//        // Calculate equivalent relative orientations, measured from particle 1
//        relOrientations[0] = q2 * q1.conj();
//        relOrientations[1] = (q2*symmetryQuaternions[0]) * q1.conj();
//        relOrientations[2] = q2 * (q1*symmetryQuaternions[0]).conj();
//        relOrientations[3] = (q2*symmetryQuaternions[0]) * (q1*symmetryQuaternions[0]).conj();
//        // Calculate equivalent relative orientations, measured from particle 2
//        relOrientations[4] = relOrientations[0].conj();
//        relOrientations[5] = relOrientations[1].conj();
//        relOrientations[6] = relOrientations[2].conj();
//        relOrientations[7] = relOrientations[3].conj();
//
//        // Looping over all equivalent relative orientations, determines if it is in state A (1) or B (2)
//        for (auto &relOrient : relOrientations) {
//            if (msmrdtools::quaternionAngleDistance(relOrient, rotMetastableStates[0]) < 2*M_PI*tolerance) {
//                return 1;
//            }
//            if (msmrdtools::quaternionAngleDistance(relOrient, rotMetastableStates[1]) < 2*M_PI*tolerance) {
//                return 2;
//            }
//        }
//        return -1; //prevsample;
//    };

//    void patchyDimer::sampleDiscreteTrajectory(double time, std::vector<particle> &particleList) {
//        std::vector<int> sample;
//        vec3<double> relativePosition; // Measured from i to j
//        quaternion<double> orientation1;
//        quaternion<double> orientation2;
//        //quaternion<double> relativeOrientation; // Measured from i to j
//        // Rotated relative position vectors w/respect to each particle
//        vec3<double> rotatedRij;
//        vec3<double> rotatedRji;
//        // Section numbers corresponding to the relative orientations in the sphere Partition
//        int secNum1;
//        int secNum2;
//
//        // Extract discrete trajectory from 2 particle simulation
//        int i = 0; // index of particle 1
//        int j = 1; // index of particle 2
//        sample = std::vector<int>{0};
//        // Extract info from two rlevant particles in list
//        orientation1 = particleList[i].orientation;
//        orientation2 = particleList[j].orientation;
//        // Calculate relative distance. If box periodic boundary, take that into account.
//        if (boundaryActive and domainBoundary->getBoundaryType() == "periodic") {
//            auto boxsize = domainBoundary->boxsize;
//            relativePosition = msmrdtools::distancePeriodicBox(particleList[j].position, particleList[i].position, boxsize);
//        } else {
//            relativePosition = particleList[j].position - particleList[i].position;
//        }
//        // Rotate the relative position vector to the frame of reference fixed in each of the two particles.
//        rotatedRij = msmrdtools::rotateVec(relativePosition, orientation1);
//        rotatedRji = msmrdtools::rotateVec(-1*relativePosition, orientation2);
//        // Evaluate if the two particles are in a bound state
//        if (relativePosition.norm() < 1.2) {
//            sample = std::vector<int>{ getBoundState(orientation1, orientation2) };
//        } else if (relativePosition.norm() < 2.0) {
//            // Get corresponding section numbers from spherical partition to classify its state
//            secNum1 = spherePart->getSectionNumber(rotatedRij);
//            secNum2 = spherePart->getSectionNumber(rotatedRji);
//            sample  = std::vector<int>{ std::min(secNum1,secNum2)*10 + std::max(secNum1,secNum2) };
//        }
//        prevsample = sample[0];
//        discreteTrajectoryData.push_back(sample);
//
//    };

//    /* Given two orientations, return either bound state A (1) or bound state B (2).
//     * Note it can be modified to take into account symmetries, using the symmetryQuaternions. */
//    int patchyDimer::getBoundState(quaternion<double> q1, quaternion<double> q2) {
//        std::array<quaternion<double>,8> relOrientations;
//        // Calculate equivalent relative orientations, measured from particle 1
//        relOrientations[0] = q2 * q1.conj();
//        relOrientations[1] = (q2*symmetryQuaternions[0]) * q1.conj();
//        relOrientations[2] = q2 * (q1*symmetryQuaternions[0]).conj();
//        relOrientations[3] = (q2*symmetryQuaternions[0]) * (q1*symmetryQuaternions[0]).conj();
//        // Calculate equivalent relative orientations, measured from particle 2
//        relOrientations[4] = relOrientations[0].conj();
//        relOrientations[5] = relOrientations[1].conj();
//        relOrientations[6] = relOrientations[2].conj();
//        relOrientations[7] = relOrientations[3].conj();
//
//        // Looping over all equivalent relative orientations, determines if it is in state A (1) or B (2)
//        for (auto &relOrient : relOrientations) {
//            if (msmrdtools::quaternionAngleDistance(relOrient, rotMetastableStates[0]) < 2*M_PI*tolerance) {
//                return 1;
//            }
//            if (msmrdtools::quaternionAngleDistance(relOrient, rotMetastableStates[1]) < 2*M_PI*tolerance) {
//                return 2;
//            }
//        }
//        return prevsample;
//    };

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