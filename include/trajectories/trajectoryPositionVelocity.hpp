//
// Created by maojrs on 10/5/21.
//
#pragma once
#include "trajectory.hpp"

namespace msmrd{
    /**
     * Class to store position and velocity trajectories
     */
    class trajectoryPositionVelocity: public trajectory {
    public:

        trajectoryPositionVelocity(unsigned long Nparticles, int bufferSize);

        void sample(double time, std::vector<particle> &particleList) override;

        void sampleRelative(double time, std::vector<particle> &particleList) override;

        // Empty function to be overwritten by child classes if neccesary
        void sampleDiscreteTrajectory(double time, std::vector<particle> &particleList) override {};

    };

    class trajectoryPositionVelocityType: public trajectoryPositionVelocity {
    public:
        using trajectoryPositionVelocity::trajectoryPositionVelocity;

        void sample(double time, std::vector<particle> &particleList) override;

        // Empty function to be overwritten by child classes if neccesary
        void sampleRelative(double time, std::vector<particle> &particleList) override {};

    };

    /*
     * It just samples the positions and velocities of distinguished particle(s) of a preselected
     * type(s). Useful for some applications, where one only needs data of
     * distinguished particle(s)
     */
    class trajectoryPositionVelocityDistinguished: public trajectoryPositionVelocityType {
    protected:
        std::vector<int> distinguishedTypes;
        /**
         * @distinguishedTypes vector of types that correspond to distinguished particle. Must match that of integrator.
         */
    public:
        trajectoryPositionVelocityDistinguished(unsigned long Nparticles, int bufferSize, std::vector<int> distinguishedTypes);

        void sample(double time, std::vector<particle> &particleList) override;

    };

}
