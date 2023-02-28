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

        void integrateO(std::vector<particle> &parts, double timestep) override;

    public:
        integratorMoriZwanzig(double dt, long seed, std::string particlesbodytype, double Gamma);

        //void integrate(std::vector<particle> &parts) override;

        void integrateOneTimestep(std::vector<particle> &parts, double timestep) override;

        void setDistinguishedTypes(std::vector<int> newDistinguishedTypes) {
            distinguishedTypes = newDistinguishedTypes; }
    };

    /*
    * Same as Mori-Zwanzig, but assumes deterministic system (zero noise).
    */
    class integratorMoriZwanzigDeterministic : public integratorMoriZwanzig {
    protected:
        void loadAuxiliaryValues(std::vector<particle> &parts, std::vector<vec3<double>> pairsForces) override;

        void integrateO(std::vector<particle> &parts, double timestep) override;

    public:
        using integratorMoriZwanzig::integratorMoriZwanzig;

    };

    /*
     * Same as integratorMoriZwanzig, but it splits the value originally store in raux to a part stored in raux and
     * another part stored in raux2.
     */
    class integratorMoriZwanzig2 : public integratorMoriZwanzig {
    protected:

        void loadAuxiliaryValues(std::vector<particle> &parts, std::vector<vec3<double>> pairsForces) override;

    public:
        using integratorMoriZwanzig::integratorMoriZwanzig;
    };

    /*
 * Same as integratorMoriZwanzig, but it splits the value originally store in raux to a part stored in raux and
 * another part stored in raux2.
 */
    class integratorMoriZwanzigConstrained1D : public integratorMoriZwanzig {
    protected:

        int index1D = 0;

        void loadAuxiliaryValues(std::vector<particle> &parts, std::vector<vec3<double>> pairsForces) override;

    public:
        using integratorMoriZwanzig::integratorMoriZwanzig;

        void integrateOneTimestep(std::vector<particle> &parts, double timestep) override;

        void updatePositionOrientation1D(std::vector<particle> &parts);

        void updateVelocities1D(std::vector<particle> &parts);


    };
}