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



    /*
     * Declaration of potential function between patchy particles with an additional explicit angular
     * dependence in terms of quaternions. Used for more precise control of angular potential over
     * orientation dynamics. (Uses the three degress of freedom of rotations to calculate torque from
     * potential, unlike its parent class uses a unit vector to describe orientation.) Used as base model
     * for pentamer formation (implemented with only one stable angular configuration).
     */
    class patchyParticleAngular2 : public patchyParticleAngular {
    public:
        // Inherit parent class constructor
        using patchyParticleAngular::patchyParticleAngular;

        double evaluate(particle &part1, particle &part2) override;

        std::array<vec3<double>, 4>
        forceTorque(particle &part1, particle &part2) override;


        // Additional auxiliary functions

        std::array<vec3<double>, 4> normalForceTorque(particle &part1, particle &part2);

        std::tuple<vec3<double>, vec3<double>, vec3<double>, vec3<double>> forceTorquePatchesSelective(
                particle &part1, particle &part2, const vec3<double> pos1virtual);

        bool isPatchBindingActive(particle &part1, particle &part2, int indexPatch1, int indexPatch2);

        std::tuple<vec3<double>, vec3<double>> calculatePlanes(const particle &part1, const particle &part2);

    };

}
