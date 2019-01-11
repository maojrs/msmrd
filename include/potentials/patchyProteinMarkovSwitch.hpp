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

        void checkMSM(particle &part1, particle &part2);
    public:
        // Inherit parent class constructor
        using patchyProtein::patchyProtein;

        double evaluate(particle &part1, particle &part2) override;

        std::array<vec3<double>, 4> forceTorque(particle &part1, particle &part2) override;

    };


}

