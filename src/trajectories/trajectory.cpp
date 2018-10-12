//
// Created by dibakma on 18.09.18.
//

#include <iostream>
#include "trajectories/trajectory.hpp"
#include "particle.hpp"


// Needed to write to HDF5 files.
using namespace H5;

namespace msmrd {

    /**
     * Implementation of abstract parent trajectory class
     */
    trajectory::trajectory(int Nparticles, int bufferSize): Nparticles(Nparticles), bufferSize(bufferSize){};

    /**
     * Implementation of trajectory class to store full position only trajectories
     */
    trajectoryPosition::trajectoryPosition(int Nparticles, int bufferSize) : trajectory(Nparticles, bufferSize){
        data.resize(0);
        data.reserve(Nparticles*bufferSize);
    };

    // Sample from list of particles and store in data
    void trajectoryPosition::sample(double time, std::vector<particle> &particleList) {
        std::array<double, 4> sample;
        for (int i = 0; i < particleList.size(); i++) {
            sample[0] = time;
            for (int k = 1; k < 4; k++) {
                sample[k] = particleList[i].position[k];
            }
            data.push_back(sample);
        }
    }

    // Sample relative positions and orientations from list of particles and store in data
    void trajectoryPosition::sampleRelative(double time, std::vector<particle> &particleList) {
        std::array<double, 4> sample;
        // Loops over all possible pairs
        for (int i = 0; i < particleList.size(); i++) {
            for (int j = i + 1; j < particleList.size(); j++) {
                sample[0] = time;
                // Relative position to particleList[j] measured from particleList[i]
                for (int k = 0; k < 3; k++) {
                    sample[k+1] = particleList[j].position[k] - particleList[i].position[k];
                }
                data.push_back(sample);
            }
        }
    };



    /**
     *  Implementation of trajectory class to store trajectories with 3D position and
     *  orientation given by a quaternion
     */
    trajectoryPositionOrientation::trajectoryPositionOrientation(int Nparticles, int bufferSize)
            : trajectory(Nparticles, bufferSize){
        data.reserve(Nparticles*bufferSize);
    };

    // Sample from list of particles and store in data
    void trajectoryPositionOrientation::sample(double time, std::vector<particle> &particleList) {
        std::array<double, 8> sample;
        for (int i = 0; i < particleList.size(); i++) {
            sample[0] = time;
            for (int k = 0; k < 3; k++) {
                sample[k+1] = particleList[i].position[k];
            }
            for (int k = 0; k < 4; k++) {
                sample[k+4] = particleList[i].orientation[k];
            }
            data.push_back(sample);
        }
    };

    // Sample relative positiona and orientation from list of particles and store in data
    void trajectoryPositionOrientation::sampleRelative(double time, std::vector<particle> &particleList) {
        std::array<double, 8> sample;
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
                data.push_back(sample);
            }
        }
    };


    void trajectoryPositionOrientation::printTime() {
        std::cerr << "Number of elements: " << data.size() << std::endl;
        for (int i=0; i<data.size(); i++) {
            std::cerr << data[i][0] << data[i][1] << data[i][2] << data[i][3] << std::endl;
        }
    };



}