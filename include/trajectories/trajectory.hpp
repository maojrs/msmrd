//
// Created by dibakma on 18.09.18.
//

#pragma once
#include <array>
#include <functional>
#include <vector>
#include "particle.hpp"

using namespace std::placeholders;

namespace msmrd {
    /**
     * Abstract base class to store full trajectories
     */
    class trajectory {
    public:
        int Nparticles;
        /**
         * @param Nparticles is the number of particles in the trajectory
         * @param data (defined in child classes) stores trajectory data (time, position, auxvariables)
         */

        trajectory(int Nparticles): Nparticles(Nparticles){};

        // Virtual function to sample from list of particles and store in data
        virtual void sample(double time, std::vector<particle> &particleList) = 0;

//        void append(sample_type sample) {
//            data.push_back(sample);
//        };
//
//        sample_type operator[](std::size_t i) const {
//            return data.at(i);
//        };
//
//        sample_type &operator[](std::size_t i) {
//            return data.at(i);
//        };
    };

    /**
     * Class to store position only trajectories
     */
    class trajectoryPosition : public trajectory {
    public:
        std::vector<std::array<double, 4>> data;

        trajectoryPosition(int Nparticles, int approx_size);

        void sample(double time, std::vector<particle> &particleList);
    };


    /**
     * Class to store trajectories with position and orientation (given by a quaternion)
     */
    class trajectoryPositionOrientation : public trajectory {
    public:
        std::vector<std::array<double, 8>> data;

        trajectoryPositionOrientation(int Nparticles, int approx_size);

        void sample(double time, std::vector<particle> &particleList);

        void printTime();
//        std::function<void(double, std::vector<particle>&)> f = std::bind(&trajectoryPositionOrientation::sample, this, _1, _2);
//        std::function<void(double, std::vector<particle>&)> get_sampler() {
//            return f;
//        }
    };


    /**
     * Class to store trajectories with relative position and relative
     * orientation (given by a quaternion) between two particles
     */
    class twoParticleRelativeTrajectory : public trajectory {
    public:
        std::vector<std::array<double, 8>> data;

        twoParticleRelativeTrajectory(int approx_size);

        void sample(double time, std::vector<particle> &particleList);
    };


} //namespace msmrd