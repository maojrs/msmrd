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
    private:
        std::array<vec3<double>, 2> refRelativePositions{};
        std::array<vec3<double>, 4> rotations{};
        std::array<quaternion<double>, 4> quatRotations{};
    public:
        // Inherit parent class constructor
        patchyParticleAngular2(double sigma, double strength, double angularStrength,
                               std::vector<std::vector<double>> patchesCoordinates);

        patchyParticleAngular2(double sigma, double strength, double angularStrength,
                               std::vector<vec3<double>> patchesCoordinates);

        using patchyParticleAngular::patchyParticleAngular;

        double evaluate(const particle &part1, const particle &part2) override;

        std::array<vec3<double>, 4>
        forceTorque(const particle &part1, const particle &part2) override;

        vec3<double> calculateQuaternionTorque(particle part1, particle part2);

        void setMetastableRegions();
    };

}

