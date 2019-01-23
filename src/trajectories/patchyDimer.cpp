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
    void patchyDimer::sampleDiscreteTrajectory(double time, std::vector<particle> &particleList) {
        int sample;
        vec3<double> relativePosition;
        quaternion<double> relativeOrientation;
        // Loops over all possible pairs (only one in this case)
        for (int i = 0; i < particleList.size(); i++) {
            for (int j = i + 1; j < particleList.size(); j++) {
                sample = 0;
                // Evaluate if the two particles are in a bound state
                relativePosition = particleList[j].position - particleList[i].position;
                relativeOrientation = particleList[j].orientation * particleList[i].orientation.conj();
                if (relativePosition.norm() <= 1.1) {
                    sample = 1;
                }
                // Relative orientation to particleList[j] measured from particleList[i]
                discreteTrajectoryData.push_back(sample);
            }
        }
    };

}