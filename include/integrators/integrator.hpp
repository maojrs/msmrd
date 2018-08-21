//
// Created by dibakma on 27.06.18.
//

#pragma once
#include <array>
#include <utility>
#include <memory>
#include "particle.hpp"
#include "randomgen.hpp"
#include "potentials/potentials.hpp"


/**
 * Abstract integrators base class declaration
 */
class integrator {
protected:
    double dt;
    long seed;
    randomgen randg;
    bool rotation;
    externalPotential *extPotential;

    /**
    * @param dt time step
    * @param seed variable for random number generation (Note seed = -1 corresponds to random device)
    * @param randg random number generator based in mt19937
    * @param rotation boolean to indicate if rotation should be integrated
    * @param extPotential potential to be used by integrator
    * @param clock keeps track of global time
    */

    /*
     * Protected abstract functions
     * Note integrateOne and integrate (public) have basically the same functionality. However, integrateOne does not
     * update the clock, so we can integrate lists and update time correctly. See implementations in src/ for details.
     */
    virtual void integrateOne(particle &part) = 0;
    virtual void translate(particle &part, double dt) = 0;
    virtual void rotate(particle &part, double dt) = 0;
public:
    double clock;

    // Base constructor
    integrator(double dt, long seed, bool rotation);;

    // Main functions definitions (=0 for abstract class)
    virtual void integrate(particle &part) = 0;
    void integrateList(std::vector<particle> &parts);
    void setExternalPotential(externalPotential *pot);
    double getClock() { return clock; }
};







