//
// Created by dibakma on 27.06.18.
//

#pragma once
#include <array>
#include <utility>
#include <memory>
#include "particle.hpp"
#include "randomgen.hpp"
#include "msm.hpp"


/**
 * Base class for integrators
 */
class integrator {
protected:
    double dt;
    long seed;
    randomgen randg;
    /**
     * @param dt time step
     * @param seed variable for random number generation;
     * @param randg random number generator based in mt19937
     */

    // Base constructor, seed=-1 corresponds to random-device seed for random number generator
    integrator(double dt, long seed): dt(dt), seed(seed) {
        randg.setSeed(seed);
    };

    // Main functions definitions (=0 for abstract class)
    virtual void integrate(particle &part) = 0;
public:
    void integrateList(std::vector<particle> &parts);
};

/**
 * Child classes of integrator
 */

// Over-damped Langevin (a.k.a. standard Brownian motion)
class odLangevin: public integrator {
public:
    odLangevin(double dt, long seed);
    void integrate(particle &part) override;

    void test(std::vector<int> &intlist);
};


// Over-damped Langevin with Markovian Switch
class odLangevinMarkovSwitch: public integrator {
public:
    std::vector<ctmsm> &msmlist;
    odLangevinMarkovSwitch(std::vector<ctmsm> &msmlist, double dt, long seed);
    void integrate(particle &part) override;
};