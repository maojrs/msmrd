//
// Created by dibakma on 18.09.18.
//

#include <iostream>
#include <fstream>
#include "trajectories/trajectory.hpp"
#include "particle.hpp"
#include "H5Cpp.h"


// Needed to write to HDF5 files.
using namespace H5;

namespace msmrd {

    /**
     * Implementation of trajectory class to store full position only trajectories
     */
    trajectoryPosition::trajectoryPosition(int Nparticles, int approx_size) : trajectory(Nparticles){
        data.reserve(Nparticles*approx_size);
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

    // Writes data into normal text file
    void trajectoryPosition::write2file(std::string filename) {
        size_t datasize = data.size();
        std::fstream outputfile;
        outputfile.open(filename, std::ios::out); // |  std::ios::binary); //vstd::ios::app |
        for(auto const& value: data)
            outputfile << value[0] << " " << value[1] << " " << value[3] << " " << value[4] << std::endl;
        //outputfile.write((char*)&data[0], datasize * sizeof(std::array<double, 4>));
        //outputfile.write((char*)&data[0], kB);
        outputfile.close();
    };

    // Writes data into HDF5 binary file
    void trajectoryPosition::write2H5file(std::string filename, std::string datasetname) {
        const H5std_string FILE_NAME = filename + ".h5";
        const H5std_string	DATASET_NAME = datasetname;
        auto writetype = H5F_ACC_RDWR; // Appends data to existing file

        // Opens a new/existing file and dataset.
        if (firstrun) {
            // CHange write type to overwrite previous file/create new file
            writetype = H5F_ACC_TRUNC;
            firstrun = false;
        }
        H5File file(FILE_NAME, writetype);

        // Sets shape of data into dataspace
        hsize_t dims[2];               // dataset dimensions
        dims[0] = data.size();
        dims[1] = 4;
        DataSpace dataspace(2, dims);

        // Creates dataset and write data into it
        DataSet dataset = file.createDataSet(DATASET_NAME, PredType::NATIVE_DOUBLE, dataspace);
        for (auto row = data.begin(); row != data.end(); row++) {
            dataset.write(&row, H5::PredType::NATIVE_DOUBLE);
        }

    };


    /**
     *  Implementation of trajectory class to store trajectories with 3D position and
     *  orientation given by a quaternion
     */
    trajectoryPositionOrientation::trajectoryPositionOrientation(int Nparticles, int approx_size)
            : trajectory(Nparticles){
        data.reserve(Nparticles*approx_size);
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




} //namespace msmrd