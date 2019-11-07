//
// Created by maojrs on 7/25/18.
//
#include <math.h>
#include "randomgen.hpp"
#include "vec3.hpp"

namespace msmrd {
    /**
     *  Implementation of random number generator class functions (constructor in header)
     */

    // Sets seed for random generator mt19937_64
    void randomgen::setSeed(long newseed) {
        if (newseed == -1) {
            mt_rand = std::mt19937_64(rd()); //random device seed
        } else if (newseed < -1) {
            mt_rand = std::mt19937_64(-newseed*rd()); //random device and time used as random seed
        } else {
            mt_rand = std::mt19937_64(newseed); //fixed seed
        }
    };

    // Returns random number between rmin and rmax sampled uniformly
    double randomgen::uniformRange(double rmin, double rmax) {
        std::uniform_real_distribution<double> uniform(rmin, rmax);
        return uniform(mt_rand);
    };

    int randomgen::uniformInteger(int imin, int imax) {
        auto value = uniformRange(1.0*imin, 1.0*imax + 1);
        return (int) std::floor(value);
    };


    // Returns random number sampled from normal distribution with given mean and standard deviation
    double randomgen::normal(double mean, double stddev) {
        std::normal_distribution<double> normaldist(mean, stddev);
        return normaldist(mt_rand);
    };

    // Returns random 3D vector sampled from 3D normal distribution with given mean and standard deviation
    vec3<double> randomgen::normal3D(double mean, double stddev) {
        vec3<double> randvec;
        std::normal_distribution<double> normaldist(mean, stddev);
        randvec[0] = normaldist(mt_rand);
        randvec[1] = normaldist(mt_rand);
        randvec[2] = normaldist(mt_rand);
        return randvec;
    };

    // Samples random 3D vector inside sphere uniformly
    vec3<double> randomgen::uniformSphere(double maxrad) {
        double rr = uniformRange(0, 1);
        double th = 2.0 * uniformRange(0, 1) - 1.0;
        double randph = 2 * M_PI * uniformRange(0, 1);
        rr = maxrad * std::pow(rr, 1. / 3.);
        th = std::acos(th);
        vec3<double> result;
        result[0] = rr * std::sin(th) * std::cos(randph);
        result[1] = rr * std::sin(th) * std::sin(randph);
        result[2] = rr * std::cos(th);
        return result;
    };

    // Samples random 3D vector inside spherical shell uniformly
    vec3<double> randomgen::uniformShell(double minrad, double maxrad) {
        double rr = -1;
        double th;
        double randph;
        while (rr < minrad || rr > maxrad) {
            rr = uniformRange(0, 1);
            rr = maxrad * std::pow(rr, 1. / 3.);
        }
        th = 2.0 * uniformRange(0, 1) - 1.0;
        th = std::acos(th);
        randph = 2 * M_PI * uniformRange(0, 1);
        vec3<double> result;
        result[0] = rr * std::sin(th) * std::cos(randph);
        result[1] = rr * std::sin(th) * std::sin(randph);
        result[2] = rr * std::cos(th);
        return result;
    }

    /* Samples 3D unit vector in sphere surface delimited by an interval in the polar angle (phiInterval)
     * and an interval in the azimuthal angle (thetaInterval) */
    vec3<double> randomgen::uniformSphereSection (std::array<double,2> polarInterval,
                                                  std::array<double,2> azimuthInterval) {
        double limTheta1 = (std::cos(polarInterval[0]) + 1.0)/2.0;
        double limTheta2 = (std::cos(polarInterval[1]) + 1.0)/2.0;
        double randomTheta = uniformRange(limTheta1, limTheta2);
        double th = std::acos(2.0 * randomTheta - 1.0);
        double randph = uniformRange(azimuthInterval[0], azimuthInterval[1]);
        vec3<double> result;
        result[0] = std::sin(th) * std::cos(randph);
        result[1] = std::sin(th) * std::sin(randph);
        result[2] = std::cos(th);
        return result;
    }

    // Samples 3D vector in a spherical shell delimited by rInterval, phiInterval and thetaInterval.
    vec3<double> randomgen::uniformShellSection (std::array<double,2> rInterval, std::array<double,2> polarInterval,
                                                 std::array<double,2> azimuthInterval) {
        // Calculate radial value (alternative to the rejection sampling used in uniformShell)
        double limR1 = std::pow(rInterval[0]/rInterval[1], 3);
        double limR2 = 1.0;
        double randomR = uniformRange(limR1, limR2);
        double rr = rInterval[1] * std::pow(randomR, 1. / 3.);
        // Calculate polar value
        double limTheta1 = (std::cos(polarInterval[0]) + 1.0)/2.0;
        double limTheta2 = (std::cos(polarInterval[1]) + 1.0)/2.0;
        double randomTheta = uniformRange(limTheta1, limTheta2);
        double th = std::acos(2.0 * randomTheta - 1.0);
        // Calculate azimuthal value
        double randph = uniformRange(azimuthInterval[0], azimuthInterval[1]);
        vec3<double> result;
        result[0] = rr * std::sin(th) * std::cos(randph);
        result[1] = rr * std::sin(th) * std::sin(randph);
        result[2] = rr * std::cos(th);
        return result;

    }

}