//
// Created by maojrs on 4/19/21.
//

#pragma once
#include "integrators/overdampedLangevin.hpp"
#include "particle.hpp"
#include "potentials/potentials.hpp"

namespace msmrd {
    /**
     * Over-damped Langevin integrator declaration. Note the initial velocity will be zero for all particles. Note
     * the integration uses the BAOAB integration scheme found in Leimkuhler's book.
     */
    class langevin : public integrator {
    protected:
        bool firstRun = true;
        /**
         * @param keeps track of first run, where double calculation of force field is required.
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

        langevin(double dt, long seed, std::string particlesbodytype);

        langevin(double dt, long seed, std::string particlesbodytype, std::string integratorScheme);

        void integrate(std::vector<particle> &parts) override;

    };

}
