//
// Created by dibakma on 17.08.18.
//

#pragma once
#include "potentials.hpp"

namespace msmrd {
    /*
     * Declaration of the Gay Berne potential which is an anisotropic Lennard Jones potential
     * Check http://www.sklogwiki.org/SklogWiki/index.php/Gay-Berne_model for details.
     * The template <vec3<double>, vec3<double>> indicates the pair-potential depends on
     * the orientations of the two particles, each described by a vector (rod-like particles).
     */
    class gayBerne : public pairPotential<vec3<double>, vec3<double>> {
    private:
        double chi;
        double chip;
        double eps0;
        double sig0;
    public:
        double a;
        double d;

        /**
         * @param a steric anisotropy of the particles: Length to breadth ratio
         * @param d ratio of binding energy in the side-by-side to end-to-end configuration
         * @param eps0 depth of the potential well
         * @param sig0 length of the short axis of the ellipsoid
         */
        gayBerne(double a, double d, double eps0, double sig0);

        double evaluate(vec3<double> pos1, vec3<double> pos2, vec3<double> u1, vec3<double> u2) override;

        std::array<vec3<double>, 2>
        forceTorque(vec3<double> pos1, vec3<double> pos2, vec3<double> u1, vec3<double> u2) override;
    };

}