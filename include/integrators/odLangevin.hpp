//
// Created by maojrs on 8/16/18.
//

#pragma once
#include "integrators/integrator.hpp"
#include "particle.hpp"
#include "potentials/potential.hpp"

/**
 * Over-damped Langevin integrator definition (a.k.a. standard Brownian motion)
 */

class odLangevin: public integrator {
protected:
    void integrateOne(particle &part) override;
    void translate(particle &part, double dt) override;
    void rotate(particle &part, double dt) override;
public:
    potential *potential0;

    odLangevin(double dt, long seed, bool rotation);

    void setPotential(potential *pot);
    void integrate(particle &part) override;
};