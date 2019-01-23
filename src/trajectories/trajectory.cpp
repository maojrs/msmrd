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
     * @param Nparticles number of particles in the simulation that need to be saved in the trajectory
     * @param bufferSize buffer size (per particle) for the trajectoryData array. If not using H5, it will be
     * used as first estimate to initialize the trajectoryData array (recommended == timeIterations/stride). If
     * using H5, it will determine the size stored in memory before flushing data into file and emptyinf buffer.
     */
    trajectory::trajectory(unsigned long Nparticles, int bufferSize): Nparticles(Nparticles), bufferSize(bufferSize){};

    // Implementation of write2file function: writes data into normal text file
    void trajectory::write2file(std::string filename, std::vector<std::vector<double>> localdata) {
        std::ofstream outputfile(filename + ".txt");
        std::ostream_iterator<double> output_iterator(outputfile, " ");

        for (auto const &value: localdata) {
            std::copy(value.begin(), value.end(), output_iterator);
            outputfile << std::endl;
        }
        outputfile.close();
    };

    /**
     * Implementation of trajectory class to store full position only trajectories
     */
    trajectoryPosition::trajectoryPosition(unsigned long Nparticles, int bufferSize) : trajectory(Nparticles, bufferSize){
        trajectoryData.resize(0);
        trajectoryData.reserve(Nparticles*bufferSize);
    };

    // Sample from list of particles and store in trajectoryData
    void trajectoryPosition::sample(double time, std::vector<particle> &particleList) {
        std::vector<double> sample(4);
        for (int i = 0; i < particleList.size(); i++) {
            sample[0] = time;
            for (int k = 1; k < 4; k++) {
                sample[k] = particleList[i].position[k];
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