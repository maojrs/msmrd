//
// Created by dibakma on 18.09.18.
//

#include <iostream>
#include "trajectory.hpp"
#include "particle.hpp"

namespace msmrd {
    trajectoryPositionOrientation::trajectoryPositionOrientation(int Nparticles, int approx_size) : trajectory(Nparticles){
        data.reserve(Nparticles*approx_size);
    };

    void trajectoryPositionOrientation::sample(double time, std::vector<msmrd::particle> &particleList) {
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

//    void trajectoryPosition::sample(double time, std::vector<msmrd::particle> &particleList) {
//        std::array<double, 4> sample;
//        for (int i = 0; i < particleList.size(); i++) {
//            sample[0] = time;
//            for (int j = 1; j < 4; j++) {
//                sample[j] = particleList[i].position[j];
//            }
//            data.push_back(sample);
//        }
//    }
} //namespace msmrd