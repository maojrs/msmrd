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

    /**
    * @param dt time step
    * @param seed variable for random number generation (Note seed = -1 corresponds to random device)
    * @param randg random number generator based in mt19937
    * @param rotation boolean to indicate if rotation should be integrated
    * @param externalPot external potential to be used by integrator
    * @param pairPot pair potential between two particles to be used by integrator
    * @param rodPairPot potential between two rod-like particles to be used by integrator
    * @param clock keeps track of global time
    * Note all potentials default to zero and not every integrator will make use of all this potentials
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
    externalPotential* externalPot;
    pairPotential* pairPot;
    rodPairPotential* rodPairPot;
    double clock;

    integrator(double dt, long seed, bool rotation);

    // Main functions definitions (=0 for abstract class)
    virtual void integrate(particle &part) = 0;
    void integrateList(std::vector<particle> &parts);
    double getClock() { return clock; }
    // Potential related functions
    void setExternalPotential(externalPotential *pot);
    void setPairPotential(pairPotential *pot);
    void setRodPairPotential(rodPairPotential *pot);
//    double evalExternalPotential(std::vector<double> pos);
//    double evalPairPotential(std::vector<double> pos1, std::vector<double> pos2);
//    double evalRodPairPotential(std::vector<double> pos1,
//                                std::vector<double> pos2,
//                                std::vector<double> u1,
//                                std::vector<double> u2);
//    std::vector<double> evalExternalForce(std::vector<double> pos);
//    std::vector<double> evalPairForce(std::vector<double> pos1, std::vector<double> pos2);
//    std::vector<std::vector<double>> evalRodPairForce(std::vector<double> pos1,
//                                                      std::vector<double> pos2,
//                                                      std::vector<double> u1,
//                                                      std::vector<double> u2);
};







