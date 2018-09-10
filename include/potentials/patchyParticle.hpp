//
// Created by maojrs on 9/10/18.
//

#pragma once
#include "potentials.hpp"

namespace msmrd {
    /*
     * Declaration of potential function between patchy particles. The template
     * <quaternion<double>,quaternion<double>> indicates the pair potential depends on the
     * orientation of both particles, each with all its three degrees of freedom.
     */
    class patchyParticle : public pairPotential<quaternion<double>, quaternion<double>> {
    private:
        std::vector<vec3<double>> patchesCoordinates;
        double sigma;

        double quadraticPotential(double r, double rstar, double eps, double a);
        double derivativeQuadraticPotential(double r, double rstar, double eps, double a);

    public:
        /*
         * @param sigma diameter of sphere at which patches are placed. Corresponds to the diameter of the particle.
         * @param patchesCoordintates list of patches coordinates in a sphere of unit radius
         */
        patchyParticle(double sigma, std::vector<vec3<double>> patchesCoordinates);
        patchyParticle(double sigma, std::vector<std::vector<double>> patchesCoordinates);


        double evaluate(vec3<double> pos1, vec3<double> pos2, quaternion<double> u1, quaternion<double> u2) override;

        std::array<vec3<double>, 2>
        forceTorque(vec3<double> pos1, vec3<double> pos2, quaternion<double> u1, quaternion<double> u2) override;

    };
}