//
// Created by maojrs on 1/14/19.
//

#pragma once
#include "potentials/patchyParticle.hpp"

namespace msmrd {
    /*
     * Declaration of potential function between patchy particles with an additional explicit angular
     * dependence.
     */
    class patchyParticleAngular : public patchyParticle {
    public:
        // Inherit parent class constructor
        using patchyParticle::patchyParticle;

        double evaluate(const particle &part1, const particle &part2) override;

        std::array<vec3<double>, 4>
        forceTorque(const particle &part1, const particle &part2) override;

    };

}
