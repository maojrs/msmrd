//
// Created by maojrs on 8/16/18.
//

#pragma once
#include "integrators/integrator.hpp"
#include "particle.hpp"
#include "potentials/potentials.hpp"

namespace msmrd {
    /**
     * Over-damped Langevin integrator declaration (a.k.a. standard Brownian motion)
     */
    class overdampedLangevin : public integrator {
    protected:
        void integrateOne(int partIndex, std::vector<particle> &parts, double timestep) override;

        void translate(particle &part, vec3<double> force, double dt) override;

        void rotate(particle &part, vec3<double> torque, double dt) override;

    public:
        overdampedLangevin(double dt, long seed, bool rotation);
    };

}