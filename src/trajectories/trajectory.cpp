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
        data.resize(0);
        data.reserve(Nparticles*bufferSize);
    };

    // Sample from list of particles and store in data
    void trajectoryPosition::sample(double time, std::vector<particle> &particleList) {
        std::vector<double> sample(4);
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
        std::vector<double> sample(4);
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
    trajectoryPositionOrientation::trajectoryPositionOrientation(unsigned long Nparticles, int bufferSize)
            : trajectory(Nparticles, bufferSize){
        data.reserve(Nparticles*bufferSize);
    };

    // Sample from list of particles and store in data
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
            data.push_back(sample);
        }
    };

    // Sample relative positiona and orientation from list of particles and store in data
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