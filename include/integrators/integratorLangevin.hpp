//
// Created by maojrs on 4/19/21.
//

#pragma once
#include "integrators/overdampedLangevin.hpp"
#include "particle.hpp"
#include "potentials/potentials.hpp"

namespace msmrd {
    /**
     * Langevin integrator declaration. The initial velocity will be zero for all particles, unless
     * explicitly given. The parent class is a virtual class. The child classes implement different algorithms
     * , e.g. the BAOAB or ABOBA schemes found in Leimkuhler's book. The different parts of some of the splitting
     * algorithms (A, B and O) are implemented in the parent class, as they are likely used by the child classes.
     * Note the integratorLangevin type integratos don't use the diffusion coefficient of the particle but the friction
     * coefficient.
     */
    class integratorLangevin : public integrator {
    protected:
        bool firstRun = true;
        double Gamma;
        /**
         * @param firstRun keeps track of first run, where double calculation of force field is required.
         * @param Gamma value of friction Coefficient for Langevin simulation, assumes the same
         * value for all particles. It must have units of mass/time
         */

        virtual void integrateOneTimestep(std::vector<particle> &parts, double timestep) = 0;

        virtual void integrateB(std::vector<particle> &parts, double timestep);

        virtual void integrateA(std::vector<particle> &parts, double timestep);

        virtual void integrateO(std::vector<particle> &parts, double timestep);

        template< typename PARTICLE >
        void updateVelocities(std::vector<PARTICLE> &parts);

        // Override empty functions (need to change integrate implementation approach for Langevin integrators)

        void integrateOne(int partIndex, std::vector<particle> &parts, double timestep) override {};

        void translate(particle &part, vec3<double> force, double dt) override {};

        void rotate(particle &part, vec3<double> torque, double dt) override {};

    public:
        integratorLangevin(double dt, long seed, std::string particlesbodytype, double Gamma);

        void integrate(std::vector<particle> &parts) override;

    };

    /* Update velocities (sets calculated next velocity
     * calculated by integrator and boundary as current velocity). Only
     * updated if particle is active and if velocityIntegration is active. */
    template <typename PARTICLE>
    void integratorLangevin::updateVelocities(std::vector<PARTICLE> &parts) {
        for (int i = 0; i < parts.size(); i++) {
            if (parts[i].isActive()) {
                parts[i].updateVelocity();
            }
        }
    }


    /*
     * Declaration of different Langevin integrators using different schemes, each in a different class
     */

    class langevinSemiImplicitEuler : public integratorLangevin {
    public:
        using integratorLangevin::integratorLangevin;

        void integrateOneTimestep(std::vector<particle> &parts, double timestep) override;

    };

    class langevinBAOAB : public integratorLangevin {
    public:
        using integratorLangevin::integratorLangevin;

        void integrateOneTimestep(std::vector<particle> &parts, double timestep) override;

    };

    class langevinABOBA : public integratorLangevin {
    public:
        using integratorLangevin::integratorLangevin;

        void integrateOneTimestep(std::vector<particle> &parts, double timestep) override;

    };

}
