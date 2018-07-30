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
    virtual void translate(particle &part, double dt) = 0;
    virtual void rotate(particle &part, double dt) = 0;
public:
    void integrateList(std::vector<particle> &parts);
};

/**
 * Child classes of integrator
 */

// Over-damped Langevin (a.k.a. standard Brownian motion)
class odLangevin: public integrator {
protected:
    void translate(particle &part, double dt) override;
    void rotate(particle &part, double dt) override;
public:
    odLangevin(double dt, long seed);
    void integrate(particle &part) override;
    double test(particle &part);
};


// Over-damped Langevin with Markovian Switch (TMSM can be an msm or a ctmsm)
template<typename TMSM>
class odLangevinMarkovSwitch: public odLangevin {
private:
    std::string msmtype;
public:
    TMSM tmsm;
    odLangevinMarkovSwitch(msm &tmsm, double dt, long seed);
    odLangevinMarkovSwitch(ctmsm &tmsm, double dt, long seed);
    void integrate(particle &part) override;
};