//
// Created by maojrs on 4/19/21.
//

#pragma once
#include "integrators/overdampedLangevin.hpp"
#include "particle.hpp"
#include "potentials/potentials.hpp"

namespace msmrd {
    /**
     * Over-damped Langevin integrator declaration. Note the initial velocity will be zero for all particles.
     */
    class langevin : public integrator {
    protected:
        void integrateOne(int partIndex, std::vector<particle> &parts, double timestep) override;

        void translate(particle &part, vec3<double> force, double dt) override;

        void rotate(particle &part, vec3<double> torque, double dt) override {};

    public:
        langevin(double dt, long seed, std::string particlesbodytype);
    };

}
