//
// Created by maojrs on 2/6/19.
//

#include "integrators/msmrdIntegrator.hpp"


namespace msmrd {

    /**
     * Implementation of base class for MSM/RD integration... to be filled
     */
    msmrdIntegrator::msmrdIntegrator(double dt, long seed, std::string particlesbodytype) :
            overdampedLangevin(dt, seed, particlesbodytype) {

    }

    void msmrdIntegrator::integrate(std::vector<particle> &parts) {}

    void msmrdIntegrator::integrateOne(int partIndex, std::vector<particle> &parts, double timestep) {};


}