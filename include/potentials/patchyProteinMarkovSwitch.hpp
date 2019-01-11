//
// Created by maojrs on 1/11/19.
//

#pragma once
#include "particle.hpp"
#include "potentials/patchyProtein.hpp"

namespace msmrd{

    /*
     * Declaration of potential function for patchy protein model with Markov switch. This class is a child from the
     * general patchy protein class that accepts quaternion-based orientation and particle types
     * rigidbody. Some of its functionality is based on the patchy particles class. Note the
     * potential will depend on the position of both particles positions, their orientation (quaternion<double>,
     * quaternion<double), and their state.
     */
    class patchyProteinMarkovSwitch : public patchyProtein {
    private:

        void setPotentialParameters() override;

        void enableDisableMSM(particleMS &part1, particleMS &part2);

    public:
        // Inherit parent class constructor
        using patchyProtein::patchyProtein;

        double evaluate(particleMS &part1, particleMS &part2);

        std::array<vec3<double>, 4> forceTorque(particleMS &part1, particleMS &part2);

        std::vector<std::vector<double>> forceTorquePyBind(particleMS &part1, particleMS &part2);

    };


}

