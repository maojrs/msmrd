//
// Created by dibakma on 18.09.18.
//

#include <iostream>
#include "trajectories/trajectory.hpp"
#include "particle.hpp"

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
            for (int j = 1; j < 4; j++) {
                sample[j] = particleList[i].position[j];
            }
            data.push_back(sample);
        }
    }


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
            for (int j = 0; j < 3; j++) {
                sample[j+1] = particleList[i].position[j];
            }
            for (int k = 0; k < 4; k++) {
                sample[k+4] = particleList[i].orientation[k];
            }
            data.push_back(sample);
        }
    };

    void trajectoryPositionOrientation::printTime() {
        std::cerr << "Number of elements: " << data.size() << std::endl;
        for (int i=0; i<data.size(); i++) {
            std::cerr << data[i][0] << data[i][1] << data[i][2] << data[i][3] << std::endl;
        }
    };


    /**
     *  Implementation of trajectory class to store relative position and
     *  orientation (given by a quaternion) between two particles
     */
    twoParticleRelativeTrajectory::twoParticleRelativeTrajectory(int approx_size)
            : trajectory(2){
        data.reserve(approx_size);
    };

    // Sample from list of particles and store in data (only works for lists of two particles)
    void twoParticleRelativeTrajectory::sample(double time, std::vector<particle> &particleList) {
        std::array<double, 8> sample;
        quaternion<double> relativeOrientation;
        sample[0] = time;
        // Relative position from particleList[0]
        for (int j = 0; j < 3; j++) {
            sample[j+1] = particleList[1].position[j] - particleList[0].position[j];
        }
        // Relative orientation from particleList[0]
        relativeOrientation = particleList[1].orientation * particleList[0].orientation.conj();
        for (int k = 0; k < 4; k++) {
            sample[k+4] = relativeOrientation[k];
        }
        data.push_back(sample);
    };


} //namespace msmrd