//
// Created by maojrs on 9/4/19.
//

#pragma once
#include "potentials/patchyParticleAngular.hpp"

namespace msmrd {
    /*
     * Declaration of potential function between patchy particles with an additional explicit angular
     * dependence in terms of quaternions. Used for more precise control of angular potential over
     * orientation dynamics. (Uses the three degress of freedom of rotations to calculate torque from
     * potential, unlike its parent class uses a unit vector to describe orientation.) Used as base model
     * for pentamer formation.
     */
    class patchyParticleAngular2 : public patchyParticleAngular {
    public:
        // Inherit parent class constructor
        using patchyParticleAngular::patchyParticleAngular;

        double evaluate(particle &part1, particle &part2) override;

        std::array<vec3<double>, 4>
        forceTorque(particle &part1, particle &part2) override;


        // Additional auxiliary functions

        std::tuple<vec3<double>, vec3<double>> calculatePlanes(const particle &part1, const particle &part2);

        };

}

