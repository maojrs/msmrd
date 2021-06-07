//
// Created by maojrs on 3/4/19.
//

#include "trajectories/trajectoryPosition.hpp"

namespace msmrd {

    /**
     * Implementation of trajectory class to store full position only trajectories
     */
    trajectoryPosition::trajectoryPosition(unsigned long Nparticles, int bufferSize) : trajectory(Nparticles, bufferSize){
        trajectoryData.reserve(Nparticles*bufferSize);
    };

    // Sample from list of particles and store in trajectoryData
    void trajectoryPosition::sample(double time, std::vector<particle> &particleList) {
        std::vector<double> sample(4);
        for (int i = 0; i < particleList.size(); i++) {
            sample[0] = time;
            for (int k = 0; k < 3; k++) {
                sample[k+1] = particleList[i].position[k];
            }
            trajectoryData.push_back(sample);
        }
    }

    // Sample relative positions and orientations from list of particles and store in trajectoryData
    void trajectoryPosition::sampleRelative(double time, std::vector<particle> &particleList) {
        std::vector<double> sample(4);
        // Loops over all possible pairs
        for (int i = 0; i < particleList.size(); i++) {
            for (int j = i + 1; j < particleList.size(); j++) {
                sample[0] = time;
                // Relative position to particleList[j] measured from particleList[i]
                for (int k = 0; k < 3; k++) {
                    sample[k+1] = particleList[j].position[k] - particleList[i].position[k];
                }
                trajectoryData.push_back(sample);
            }
        }
    };


    // Sample from list of particles and store in trajectoryData (for trajectoryPositionType)
    void trajectoryPositionType::sample(double time, std::vector<particle> &particleList) {
        std::vector<double> sample(5);
        for (int i = 0; i < particleList.size(); i++) {
            sample[0] = time;
            for (int k = 0; k < 3; k++) {
                sample[k+1] = particleList[i].position[k];
            }
            sample[4] = 1.0 * particleList[i].type;
            trajectoryData.push_back(sample);
        }
    }

}
