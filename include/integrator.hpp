//
// Created by dibakma on 27.06.18.
//

#pragma once
#include <array>
#include <utility>
#include <memory>
#include "particle.hpp"
#include "randomgen.hpp"


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
    virtual void integrate(std::vector<std::shared_ptr<particle>> parts) = 0;
};

/**
 * Over-damped Langevin (a.k.a. standard Brownian motion)
 */

class odLangevin: public integrator {
public:
    odLangevin(double dt, long seed);
    void integrate(std::vector<std::shared_ptr<particle>> parts) override;
    double test();
};