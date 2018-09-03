//
// Created by maojrs on 8/16/18.
//

#pragma once
#include "integrators/odLangevin.hpp"
#include "particle.hpp"
#include "msm.hpp"

/**
 * Over-damped Langevin with Markovian Switch integrator
 * @tparam TMSM template can be an msm or a ctmsm
 */
template<typename TMSM>
class odLangevinMarkovSwitch: public odLangevin {
private:
    std::string msmtype;
    void integrateOne(particleMS &part, double timestep);
    void integrateOne(int partIndex, std::vector<particleMS> &parts, double timestep);
    void integrateOneMS(particleMS &part, double timestep);
    void integrateOneMS(int partIndex, std::vector<particleMS> &parts, double timestep);
public:
    TMSM tmsm;
    /**
     * @param msmtype string to distinguish between msm and ctmsm
     * @param tmsm object variable, can be either an msm or a ctmsm
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
    void integrateList(std::vector<particleMS> &parts, bool pairInteractions);
};