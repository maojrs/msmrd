//
// Created by maojrs on 2/6/19.
//

#pragma once
#include "integrators/integrator.hpp"
#include "overdampedLangevinMarkovSwitch.hpp"

namespace msmrd {
    /**
     * Base class for msmrd integration (coupling MSM and reaction-diffusion)
     */
    class msmrdIntegrator : public overdampedLangevin {
    protected:
        void integrateOne(int partIndex, std::vector<particle> &parts, double timestep) override;

    public:
        msmrdIntegrator(double dt, long seed, std::string particlesbodytype);

        // Redefine integrate function
        void integrate(std::vector<particle> &parts);
    };


}
