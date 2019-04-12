//
// Created by maojrs on 2/6/19.
//

#pragma once
#include "discretizations/positionOrientationPartition.hpp"
#include "integrators/integrator.hpp"
#include "markovModels/msmrdMarkovModel.hpp"
#include "overdampedLangevinMarkovSwitch.hpp"

namespace msmrd {
    using msmrdMSM = msmrd::msmrdMarkovModel;
    using fullPartition = msmrd::positionOrientationPartition;

    /**
     * Base class for msmrd integration (coupling MSM and reaction-diffusion)
     */
    class msmrdIntegrator : public overdampedLangevin {
    private:
        msmrdMSM &markovModel;
        fullPartition &positionOrientationPart;

        void integrateOne(int partIndex, std::vector<particle> &parts, double timestep) override;

        vec3<double> calculateRelativePosition(particle &p1, particle &p2);
    public:
        /**
        * @param markovModel pointer to class msmrdMarkovModel, which is the markovModel class specialized for
        * the MSM/RD scheme.
        */

        msmrdIntegrator(double dt, long seed, std::string particlesbodytype, msmrdMSM markovModel,
                        fullPartition positionOrientationPart);

        // Redefine integrate function
        void integrate(std::vector<particle> &parts, msmrdMSM &masterMSM);
    };


}
