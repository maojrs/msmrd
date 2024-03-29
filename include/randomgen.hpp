//
// Created by maojrs on 7/25/18.
//
#pragma once
#include <random>
#include "vec3.hpp"
#include <chrono>

namespace msmrd {
    /**
     * Declares class to generate good random numbers (mt_19937_64 generator)
     * and sample from common distirbutions
     */
    class randomgen {
    private:
        long seed = -1;
        std::random_device rd;
        std::mt19937_64 mt_rand;  // random number generator
    public:
        randomgen() { mt_rand = std::mt19937_64(rd()); };

        randomgen(const randomgen &) { mt_rand = std::mt19937_64(rd()); };

        void setSeed(long newseed);

        double uniformRange(double rmin, double rmax);

        int uniformInteger(int imin, int imax);

        double normal(double mean, double stddev);

        vec3<double> normal3D(double mean, double stddev);

        vec3<double> uniformSphere(double maxrad);

        vec3<double> uniformShell(double minrad, double maxrad);

        vec3<double> uniformSphereSection (std::array<double,2> polarInterval, std::array<double,2> azimuthInterval);

        vec3<double> uniformShellSection (std::array<double,2> rInterval, std::array<double,2> polarInterval,
                                          std::array<double,2> azimuthInterval);



    };

}
