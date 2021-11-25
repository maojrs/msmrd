//
// Created by maojrs on 11/25/21.
//
#pragma once
#include <memory>
#include "trajectory.hpp"
#include "potentials/potentials.hpp"

namespace msmrd {
    /**
     * Class to store energy and instant temperature of a particles in a given simulation. The potentials should be
     * set to be exactly the same ones used by the integrator.
     */
    class trajectoryEnergyTemperature : public trajectory {
    public:
        bool externalPotentialActive = false;
        bool pairPotentialActive = false;

        // External potentials pointers (externalPot)
        std::shared_ptr<externalPotential> externalPot;

        // Pair potentials pointers (pairPot)
        std::shared_ptr<pairPotential> pairPot;

        trajectoryEnergyTemperature(int bufferSize);

        void setExternalPotential(externalPotential *potential);

        void setPairPotential(pairPotential *potential);

        void sample(double time, std::vector<particle> &particleList) override;

        // Empty functions, can be overwritten by child classes if neccesary
        void sampleRelative(double time, std::vector<particle> &particleList) override {};

        void sampleDiscreteTrajectory(double time, std::vector<particle> &particleList) override {};

    };

}