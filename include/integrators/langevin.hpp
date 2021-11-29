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
     * explicitly given. The integration uses the BAOAB or ABOBA integration schemes found in Leimkuhler's book.
     * Note the langevin type integratos don't use the diffusion coefficient of the particle but the friction
     * coefficient.
     */
    class langevin : public integrator {
    protected:
        bool firstRun = true;
        double frictionCoefficient;
        /**
         * @param firstRun keeps track of first run, where double calculation of force field is required.
         * @param frictionCoefficient value of friction Coefficient for Langevin simulation, assumes the same
         * value for all particles. It must have units of mass/time
         */

        void integrateB(std::vector<particle> &parts, double timestep);

        void integrateA(std::vector<particle> &parts, double timestep);

        void integrateO(std::vector<particle> &parts, double timestep);

        void integrateBAOAB(std::vector<particle> &parts, double timestep);

        void integrateABOBA(std::vector<particle> &parts, double timestep);

        // Override empty functions (need to change integrate implementation approach for Langevin integrators)

        void integrateOne(int partIndex, std::vector<particle> &parts, double timestep) override {};

        void translate(particle &part, vec3<double> force, double dt) override {};

        void rotate(particle &part, vec3<double> torque, double dt) override {};

    public:
        std::string integratorScheme;

        langevin(double dt, long seed, std::string particlesbodytype, double frictionCoefficient);

        langevin(double dt, long seed, std::string particlesbodytype, double frictionCoefficient,
                std::string integratorScheme);

        void integrate(std::vector<particle> &parts) override;

    };

}
