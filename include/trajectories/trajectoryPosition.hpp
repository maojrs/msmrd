//
// Created by maojrs on 3/4/19.
//

#pragma once
#include "trajectory.hpp"

namespace msmrd{
    /**
     * Class to store position only trajectories
     */
    class trajectoryPosition: public trajectory {
    public:

        trajectoryPosition(unsigned long Nparticles, int bufferSize);

        void sample(double time, std::vector<particle> &particleList) override;

        void sampleRelative(double time, std::vector<particle> &particleList) override;

        // Empty function to be overwritten by child classes if neccesary
        void sampleDiscreteTrajectory(double time, std::vector<particle> &particleList) override {};

    };

    class trajectoryPositionType: public trajectoryPosition {
    public:
        using trajectoryPosition::trajectoryPosition;

        void sample(double time, std::vector<particle> &particleList) override;

        // Empty function to be overwritten by child classes if neccesary
        void sampleRelative(double time, std::vector<particle> &particleList) override {};

    };

    /*
     * It just samples the position of distinguished particle(s) of a preselected
     * type(s). Useful for some applications, where one only needs data of
     * distinguished particle(s)
     */
    class trajectoryPositionDistinguished: public trajectoryPositionType {
    protected:
        std::vector<int> distinguishedTypes;
        /**
         * @distinguishedTypes vector of types that correspond to distinguished particle. Must match that of integrator.
         */
    public:
        trajectoryPositionDistinguished(unsigned long Nparticles, int bufferSize, std::vector<int> distinguishedTypes);

        void sample(double time, std::vector<particle> &particleList) override;

    };

}