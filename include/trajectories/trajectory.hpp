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
    protected:
        const std::size_t kB = 1024;
        const std::size_t MB = 1024 * kB;
        bool firstrun = true;
    public:
        int Nparticles;
        /**
         * @param kB/MB, constant buffer sizes in bytes for writing data
         * @param Nparticles is the number of particles in the trajectory. In the case of relative
         * sampling of cooridnates, it should correspond to the number of all possible pairs of
         * particles.
         * @param data (defined in child classes to avoid template) stores
         * trajectory data (time, position, and/or other variables like orientation)
         */

        trajectory(int Nparticles): Nparticles(Nparticles){};

        // Virtual functions to sample from list of particles and store in data
        virtual void sample(double time, std::vector<particle> &particleList) = 0;

        virtual void sampleRelative(double time, std::vector<particle> &particleList) = 0;


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
    private:
        std::vector<std::array<double, 4>> data;
    public:

        trajectoryPosition(int Nparticles, int approx_size);

        void sample(double time, std::vector<particle> &particleList) override;

        void sampleRelative(double time, std::vector<particle> &particleList) override;

        void write2file(std::string filename);

        void write2H5file(std::string filename, std::string datasetname);

        std::vector<std::array<double, 4>> getData() const { return data; };

    };


    /**
     * Class to store trajectories with position and orientation (given by a quaternion)
     */
    class trajectoryPositionOrientation : public trajectory {
    private:
        std::vector<std::array<double, 8>> data;
    public:

        trajectoryPositionOrientation(int Nparticles, int approx_size);

        void sample(double time, std::vector<particle> &particleList) override;

        void sampleRelative(double time, std::vector<particle> &particleList) override;

        std::vector<std::array<double, 8>> getData() const { return data; };

        void printTime();
//        std::function<void(double, std::vector<particle>&)> f = std::bind(&trajectoryPositionOrientation::sample, this, _1, _2);
//        std::function<void(double, std::vector<particle>&)> get_sampler() {
//            return f;
//        }
    };



} //namespace msmrd