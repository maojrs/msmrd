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
    bool rotation;

    // Protected abstract functions
    virtual void translate(particle &part, double dt) = 0;
    virtual void rotate(particle &part, double dt) = 0;
public:
    double clock;
    /**
     * @param dt time step
     * @param seed variable for random number generation;
     * @param randg random number generator based in mt19937
     * @param rotation boolean to indicate if rotation should be integrated
     * @param clock keeps track of global time
     */

    // Base constructor, seed=-1 corresponds to random-device seed for random number generator
    integrator(double dt, long seed, bool rotation)
            : dt(dt), seed(seed), rotation(rotation) {
        randg.setSeed(seed);
        clock = 0;
    };

    // Main functions definitions (=0 for abstract class)
    virtual void integrate(particle &part) = 0;
    void integrateList(std::vector<particle> &parts);

    /**
    * Get and set functions. Some used by c++ and python,
    * some only to be used by pyhon with python bindings.
    **/
    double getClock() { return clock; }
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
    odLangevin(double dt, long seed, bool rotation);
    void integrate(particle &part) override;
};


// Over-damped Langevin with Markovian Switch (TMSM template can be an msm or a ctmsm)
template<typename TMSM>
class odLangevinMarkovSwitch: public odLangevin {
private:
    std::string msmtype;
public:
    TMSM tmsm;
    /**
     * @param msmtype string to distinguish between msm and ctmsm
     * @param tmsm object variable, cn be either an msm or a ctmsm
     */

    // Constructors need to be defined in headers for template w/pybind
    odLangevinMarkovSwitch(ctmsm &tmsm, double dt, long seed, bool rotation)
            : tmsm(tmsm), odLangevin(dt,seed,rotation) {
        msmtype = "continuous-time";
    };

    odLangevinMarkovSwitch(msm &tmsm, double dt, long seed, bool rotation)
            : tmsm(tmsm), odLangevin(dt,seed,rotation) {
        msmtype = "discrete-time";
    };

    void integrate(particleMS &part);
    void integrateList(std::vector<particleMS> &parts);

};
