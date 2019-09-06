//
// Created by maojrs on 1/14/19.
//

#pragma once
#include "potentials/patchyParticle.hpp"
#include "quaternion.hpp"
#include "tools.hpp"

namespace msmrd {
    /*
     * Declaration of potential function between patchy particles with an additional explicit angular
     * dependence.
     */
    class patchyParticleAngular : public patchyParticle {
    protected:
        double patchPotentialScaling = 0.5; // 0.3 original value
        double angularStrength = 2.0;
    public:
        // Inherit parent class constructor
        using patchyParticle::patchyParticle;

        // Additionals constructors in case angular strength is provided
        patchyParticleAngular(double sigma, double strength, double angularStrength,
                              std::vector<std::vector<double>> patchesCoordinates);

        patchyParticleAngular(double sigma, double strength, double angularStrength,
                              std::vector<vec3<double>> patchesCoordinates);

        double evaluate(const particle &part1, const particle &part2) override;

        std::array<vec3<double>, 4>
        forceTorque(const particle &part1, const particle &part2) override;

    };

}
