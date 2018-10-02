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
    class trajectory {
    public:
        int Nparticles;
        trajectory(int Nparticles): Nparticles(Nparticles){};

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

//    class trajectoryPosition : public trajectory<std::array<double, 4>> {
//    public:
//        void sample(double time, std::vector<particle> &particleList);
//    };
} //namespace msmrd