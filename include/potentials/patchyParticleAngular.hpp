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
        double angularStrength = 2.0;
    public:
        /**
         * @param angularStrength give the angular strength of angular dependence of torque.
         */

        // Inherit parent class constructor (in case no angular strength is provided)
        using patchyParticle::patchyParticle;

        // Additionals constructors in case angular strength is provided
        patchyParticleAngular(double sigma, double strength, double angularStrength,
                              std::vector<std::vector<double>> patchesCoordinates);

        patchyParticleAngular(double sigma, double strength, double angularStrength,
                              std::vector<vec3<double>> patchesCoordinates);

        double evaluate(particle &part1, particle &part2) override;

        std::array<vec3<double>, 4>
        forceTorque(particle &part1, particle &part2) override;

    };

}
