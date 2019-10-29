//
// Created by maojrs on 8/16/18.
//

#pragma once
#include "integrators/overdampedLangevin.hpp"
#include "particle.hpp"
#include "markovModels/discreteTimeMarkovModel.hpp"
#include "markovModels/continuousTimeMarkovModel.hpp"


namespace msmrd {
    using msm = msmrd::discreteTimeMarkovStateModel;
    using ctmsm = msmrd::continuousTimeMarkovStateModel;
    /**
     * Over-damped Langevin with Markovian Switch integrator
     * @tparam templateMSM template can be an msm or a ctmsm
     */
    template<typename templateMSM>
    class overdampedLangevinMarkovSwitch : public overdampedLangevin {
    protected:
        std::string msmtype;

        void integrateOne(int partIndex, std::vector<particle> &parts, double timestep);

        void integrateOneMS(int partIndex, std::vector<particle> &parts, double timestep);

    public:
        std::vector<templateMSM> MSMlist;
        /**
         * @param msmtype string to distinguish between msm and ctmsm
         * @param MSMlist can be either a vector of msms or of ctmsms. It is possible to provide only one msm/ctmsm (if
         * only one type of particles is being used) in constructor instead of a vector.
         */

        /* Constructors need to be defined in headers for template w/pybind, see parent
         * class overdampedLangevin for details on constructor parameters */

        overdampedLangevinMarkovSwitch(std::vector<ctmsm> &MSMlist, double dt, long seed, std::string particlesbodytype)
                : MSMlist(MSMlist), overdampedLangevin(dt, seed, particlesbodytype) {
            msmtype = "continuous-time";
        };

        overdampedLangevinMarkovSwitch(std::vector<msm> &MSMlist, double dt, long seed, std::string particlesbodytype)
                : MSMlist(MSMlist), overdampedLangevin(dt, seed, particlesbodytype) {
            msmtype = "discrete-time";
        };

        overdampedLangevinMarkovSwitch(ctmsm &MSMlist, double dt, long seed, std::string particlesbodytype)
                : MSMlist(std::vector<ctmsm>{MSMlist}), overdampedLangevin(dt, seed, particlesbodytype) {
            msmtype = "continuous-time";
        };

        overdampedLangevinMarkovSwitch(msm &MSMlist, double dt, long seed, std::string particlesbodytype)
                : MSMlist(std::vector<msm>{MSMlist}), overdampedLangevin(dt, seed, particlesbodytype) {
            msmtype = "discrete-time";
        };

        void integrate(std::vector<particle> &parts);
    };

}