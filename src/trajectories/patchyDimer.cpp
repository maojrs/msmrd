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
    patchyDimer::patchyDimer(unsigned long Nparticles, int bufferSize, int numSections) :
    trajectoryPositionOrientation(Nparticles, bufferSize) {
        spherePart = std::make_unique<spherePartition>(numSections);
        angularStates = numSections*(numSections + 1)/2;
        setMetastableRegions();
    };


    void patchyDimer::sampleDiscreteTrajectory(double time, std::vector<particle> &particleList) {
        std::vector<int> sample;
        vec3<double> relativePosition; // Measured from i to j
        quaternion<double> orientation1;
        quaternion<double> orientation2;
        //quaternion<double> relativeOrientation; // Measured from i to j
        // Rotated relative position vectors w/respect to each particle
        vec3<double> rotatedRij;
        vec3<double> rotatedRji;
        // Section numbers corresponding to the relative orientations in the sphere Partition
        int secNum1;
        int secNum2;

        // Loops over all possible pairs (only one in this case)
        for (int i = 0; i < particleList.size(); i++) {
            for (int j = i + 1; j < particleList.size(); j++) {
                sample = std::vector<int>{0};
                // Extract info from two rlevant particles in list
                orientation1 = particleList[i].orientation;
                orientation2 = particleList[j].orientation;
                relativePosition = particleList[j].position - particleList[i].position;
                // Rotate the relative position vector to the frame of reference fixed in each of the two particles.
                rotatedRij = msmrdtools::rotateVec(relativePosition, orientation1);
                rotatedRji = msmrdtools::rotateVec(-1*relativePosition, orientation2);
                // Evaluate if the two particles are in a bound state
                if (relativePosition.norm() <= 1.1) {
                    sample = std::vector<int>{ getBoundState(orientation1, orientation2) };
                } else if (relativePosition.norm() <= 2.1) {
                    // Get corresponding section numbers from spherical partition to classify its state
                    secNum1 = spherePart->getSectionNumber(rotatedRij);
                    secNum2 = spherePart->getSectionNumber(rotatedRji);
                    sample  = std::vector<int>{ std::min(secNum1,secNum2)*10 + std::max(secNum1,secNum2) };
                }
                discreteTrajectoryData.push_back(sample);
            }
        }
    };

    /* Given two orientations, return either bound state A (1) or bound state B (2).
     * Note it needs to take into account all the symmetries. */
    int patchyDimer::getBoundState(quaternion<double> q1, quaternion<double> q2) {
        std::array<quaternion<double>,8> relOrientations;
        // Calculate equivalent relative orientations, measured from particle 1
        relOrientations[0] = q2 * q1.conj();
        relOrientations[1] = (q2*symmetryQuaternion) * q1.conj();
        relOrientations[2] = q2 * (q1*symmetryQuaternion).conj();
        relOrientations[3] = (q2*symmetryQuaternion) * (q1*symmetryQuaternion).conj();
        // Calculate equivalent relative orientations, measured from particle 2
        relOrientations[4] = relOrientations[0].conj();
        relOrientations[5] = relOrientations[1].conj();
        relOrientations[6] = relOrientations[2].conj();
        relOrientations[7] = relOrientations[3].conj();

        // Looping over all equivalent relative prientations, determins if it is in state A (1) or B (2)
        for (auto &relOrient : relOrientations) {
            if (msmrdtools::quaternionDistance(relOrient, rotMetastableStates[0]) < tolerance) {
                return 1;
            }
            if (msmrdtools::quaternionDistance(relOrient, rotMetastableStates[1]) < tolerance) {
                return 2;
            }
        }
        return 0;
    };

    /*
     * Define metastable relative orientations (including symmetric quaternion)
     */
    void patchyDimer::setMetastableRegions() {
        // Define symmetry quaternion that denotes rotational symmetry of particle
        vec3<double> axisAngleRotSym = vec3<double>(M_PI, 0, 0);
        symmetryQuaternion = msmrdtools::axisangle2quaternion(axisAngleRotSym);
        // Define 2 main rotation quaternions, qmetaA and qmetaB, defining bound state A and B state.
        int numStates = 2;
        vec3<double> axisAngleRotA = vec3<double>(0, 0, M_PI - 3.0*M_PI/5);
        vec3<double> axisAngleRotB = vec3<double>(0, 0, M_PI);
        quaternion<double> qmetaA = msmrdtools::axisangle2quaternion(axisAngleRotA);
        quaternion<double> qmetaB = msmrdtools::axisangle2quaternion(axisAngleRotB);

        rotMetastableStates.resize(numStates);
        rotMetastableStates[0] = 1*qmetaA;
        rotMetastableStates[1] = 1*qmetaB;
    }

}