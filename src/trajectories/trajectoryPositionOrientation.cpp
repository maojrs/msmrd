//
// Created by maojrs on 3/4/19.
//

#include "trajectories/trajectoryPositionOrientation.hpp"

namespace msmrd {
    /**
    *  Implementation of trajectory class to store trajectories with 3D position and
    *  orientation given by a quaternion
    */
    trajectoryPositionOrientation::trajectoryPositionOrientation(unsigned long Nparticles, int bufferSize)
            : trajectory(Nparticles, bufferSize){
        trajectoryData.reserve(Nparticles*bufferSize);
    };

    // Sample from list of particles and store in trajectoryData
    void trajectoryPositionOrientation::sample(double time, std::vector<particle> &particleList) {
        std::vector<double> sample(8);
        for (int i = 0; i < particleList.size(); i++) {
            sample[0] = time;
            for (int k = 0; k < 3; k++) {
                sample[k+1] = particleList[i].position[k];
            }
            for (int k = 0; k < 4; k++) {
                sample[k+4] = particleList[i].orientation[k];
            }
            trajectoryData.push_back(sample);
        }
    };

    // Sample relative positiona and orientation from list of particles and store in trajectoryData
    void trajectoryPositionOrientation::sampleRelative(double time, std::vector<particle> &particleList) {
        std::vector<double> sample(8);
        quaternion<double> relativeOrientation;
        // Loops over all possible pairs
        for (int i = 0; i < particleList.size(); i++) {
            for (int j = i + 1; j < particleList.size(); j++) {
                sample[0] = time;
                // Relative position to particleList[j] measured from particleList[i]
                for (int k = 0; k < 3; k++) {
                    sample[k+1] = particleList[j].position[k] - particleList[i].position[k];
                }
                // Relative orientation to particleList[j] measured from particleList[i]
                relativeOrientation = particleList[j].orientation * particleList[i].orientation.conj();
                for (int k = 0; k < 4; k++) {
                    sample[k+4] = relativeOrientation[k];
                }
                trajectoryData.push_back(sample);
            }
        }
    };


    void trajectoryPositionOrientation::printTime() {
        std::cerr << "Number of elements: " << trajectoryData.size() << std::endl;
        for (int i=0; i<trajectoryData.size(); i++) {
            std::cerr << trajectoryData[i][0] << trajectoryData[i][1] << trajectoryData[i][2] << trajectoryData[i][3] << std::endl;
        }
    };

}

