//
// Created by maojrs on 6/14/21.
//

#pragma once
#include "integrators/integratorLangevin.hpp"
#include "particle.hpp"
#include "potentials/potentials.hpp"

namespace msmrd {
    /**
     * Integrator class for integratorMoriZwanzig. It is essentially the same as the integratorLangevin integrator, but
     * it stores specific values into the raux variable of the particle class for further analysis. It is also
     * only implemented for the ABOBA since it is ideal to extract the raux variables.
     *
     * It assumes the distinguished particle is the one with index0 (type 1).
     */
    class integratorMoriZwanzig : public langevinABOBA {
    protected:
        std::vector<int> distinguishedTypes{1};
        /**
         * @distinguishedTypes vector of types that correspond to distinguished particle. Must match that of trajectory.
         */
        virtual void loadAuxiliaryValues(std::vector<particle> &parts, std::vector<vec3<double>> pairsForces);

        virtual void integrateO(std::vector<particle> &parts, double timestep);

    public:
        integratorMoriZwanzig(double dt, long seed, std::string particlesbodytype, double Gamma);

        void integrate(std::vector<particle> &parts) override;

        void setDistinguishedTypes(std::vector<int> newDistinguishedTypes) {
            distinguishedTypes = newDistinguishedTypes; }



    };

    /*
     * Same as integratorMoriZwanzig, but it splits the value originally store in raux to a part stored in raux and
     * another part stored in raux2.
     */
    class integratorMoriZwanzig2 : public integratorMoriZwanzig {
    protected:
        /**
         * @distinguishedTypes vector of types that correspond to distinguished particle. Must match that of trajectory.
         */
        void loadAuxiliaryValues(std::vector<particle> &parts, std::vector<vec3<double>> pairsForces) override;

    public:
        using integratorMoriZwanzig::integratorMoriZwanzig;
    };
}