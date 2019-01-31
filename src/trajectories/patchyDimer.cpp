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
    };


    void patchyDimer::sampleDiscreteTrajectory(double time, std::vector<particle> &particleList) {
        int sample;
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
                sample = 0;
                // Extract info from two rlevant particles in list
                orientation1 = particleList[i].orientation;
                orientation2 = particleList[j].orientation;
                relativePosition = particleList[j].position - particleList[i].position;
                // Rotate the relative position vector to the frame of reference fixed in each of the two particles.
                rotatedRij = msmrdtools::rotateVec(relativePosition, orientation1);
                rotatedRji = msmrdtools::rotateVec(-1*relativePosition, orientation2);
                // Get corresponding section numbers from spherical partition to classify its state
                secNum1 = spherePart->getSectionNumber(rotatedRij);
                secNum2 = spherePart->getSectionNumber(rotatedRji);
                //relativeOrientation = orientation2 * orientation1.conj();
                // Evaluate if the two particles are in a bound state
                if (relativePosition.norm() <= 1.1) {
                    sample = 1;
                }
                // Relative orientation to particleList[j] measured from particleList[i]
                discreteTrajectoryData.push_back(sample);
            }
        }
    };


    /* MAYBE NOT USEFUL, WORK IN PROGRESS
     * Define metastable relative orientations (including symmetric representations)
     */
    void patchyDimer::setMetastableRegions() {
        int numStates = 4;
        int numSymmetries = 4;
        // Define 4 main rotation quaternions, q1, q2, q3 and q4, each defining a possible different  metstable state
        vec3<double> axisAngleRot1 = vec3<double>(0, 0, M_PI - 3.0*M_PI/5);
        vec3<double> axisAngleRot2 = vec3<double>(0, 0, -M_PI);
        vec3<double> axisAngleRot3 = vec3<double>(0, 0, 3.0*M_PI/5);
        vec3<double> axisAngleRot4 = vec3<double>(0, 0, 3.0*M_PI/5);
        quaternion<double> q1 = msmrdtools::axisangle2quaternion(axisAngleRot1);
        quaternion<double> q2;
        quaternion<double> q3;
        quaternion<double> q4;

        rotMetastableStates.resize(numStates);
        for (auto rotMstate : rotMetastableStates) {
            rotMstate.resize(numSymmetries);
            // Define equivalent rotations
            q1 = quaternion<double> ();
            q2 = quaternion<double> ();
            q3 = quaternion<double> ();
            q4 = quaternion<double> ();
            // Save equivalent rotation in corresponding list element
            rotMstate = std::vector<quaternion<double>> {q1, q2, q3, q4};
        }
    }

}