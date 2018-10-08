//
// Created by maojrs on 8/16/18.
//

#pragma once
#include "integrators/overdampedLangevin.hpp"
#include "particle.hpp"
#include "msm.hpp"

namespace msmrd {
    using msm = msmrd::discreteTimeMarkovStateModel;
    using ctmsm = msmrd::continuousTimeMarkovStateModel;
    /**
     * Over-damped Langevin with Markovian Switch integrator
     * @tparam TMSM template can be an msm or a ctmsm
     */
    template<typename TMSM>
    class overdampedLangevinMarkovSwitch : public overdampedLangevin {
    private:
        std::string msmtype;

        void integrateOne(int partIndex, std::vector<particleMS> &parts, double timestep);

        void integrateOneMS(int partIndex, std::vector<particleMS> &parts, double timestep);

    public:
        TMSM tmsm;
        /**
         * @param msmtype string to distinguish between msm and ctmsm
         * @param tmsm object variable, can be either an msm or a ctmsm
         */

        /* Constructors need to be defined in headers for template w/pybind, see parent
         * class odLangevin for details on constructor parameters */
        overdampedLangevinMarkovSwitch(ctmsm &tmsm, double dt, long seed, std::string particlesbodytype)
                : tmsm(tmsm), overdampedLangevin(dt, seed, particlesbodytype) {
            msmtype = "continuous-time";
        };

        overdampedLangevinMarkovSwitch(msm &tmsm, double dt, long seed, std::string particlesbodytype)
                : tmsm(tmsm), overdampedLangevin(dt, seed, particlesbodytype) {
            msmtype = "discrete-time";
        };

        void integrate(std::vector<particleMS> &parts);
    };

}