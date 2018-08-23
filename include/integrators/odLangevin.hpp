//
// Created by maojrs on 8/16/18.
//

#pragma once
#include "integrators/integrator.hpp"
#include "particle.hpp"
#include "potentials/potentials.hpp"

/**
 * Over-damped Langevin integrator declaration (a.k.a. standard Brownian motion)
 */
class odLangevin: public integrator {
protected:
    void integrateOne(particle &part) override;
    void translate(particle &part, vec3<double> force, double dt) override;
    void rotate(particle &part, vec3<double> torque, double dt) override;
public:
    odLangevin(double dt, long seed, bool rotation);

    void integrate(particle &part) override;
};