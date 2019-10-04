//
// Created by maojrs on 9/10/18.
//

#pragma once
#include "potentials.hpp"

namespace msmrd {
    /*
     * Declaration of potential function between patchy particles.This pair potential depends on the
     * orientation of both particles, each with all its three rotational degrees of freedom.
     */
    class patchyParticle : public pairPotential {
    protected:
        std::vector<vec3<double>> patchesCoordinates;
        double sigma = 1.0;
        double strength = 100.0;
        double patchPotentialScaling = 1.0;
        // Strength of potentials
        double epsRepulsive;
        double epsAttractive;
        double epsPatches;
        // Stiffness of potentials
        double aRepulsive;
        double aAttractive;
        double aPatches;
        // Range parameter of potentials
        double rstarRepulsive;
        double rstarAttractive;
        double rstarPatches;

        double quadraticPotential(double r, double sig, double eps, double a, double rstar);
        double derivativeQuadraticPotential(double r, double sig, double eps, double a, double rstar);
        std::tuple<vec3<double>, vec3<double>, vec3<double>, vec3<double>> forceTorquePatches(
                const particle &part1, const particle &part2, const vec3<double> pos1virtual);

    public:
        /*
         * @param patchesCoordintates list of patches coordinates in a sphere of unit radius
         * @param sigma diameter of sphere at which patches are placed. Corresponds to the diameter of the particle.
         * @param strength overall strength of potential, eps values scaled in terms of this value.
         * @param eps*** strength paramaters for isotropic attractive, isotropic repulsive and
         * anisotropic attractive patches potentials.
         * @param a*** stiffness paramaters for the same three potentials.
         * @param rstar*** range paramaters for the same three potentials.
         */
        patchyParticle() = default;
        patchyParticle(double sigma, double strength, std::vector<vec3<double>> patchesCoordinates);
        patchyParticle(double sigma, double strength, std::vector<std::vector<double>> patchesCoordinates);


        double evaluate(const particle &part1, const particle &part2) override;

        std::array<vec3<double>, 4>
        forceTorque(const particle &part1, const particle &part2) override;

    };
}