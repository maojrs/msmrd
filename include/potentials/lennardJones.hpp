//
// Created by maojrs on 4/19/21.
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
    class lennardJones : public pairPotential {
    public:
        double epsilon;
        double sigma;
        double cutOff;

        /**
         * @param epsilon depth of the potential well
         * @param sigma distance at which the potential is zero (a.k.a 'size of the particle')
         */
        lennardJones(double epsilon, double sigma);

        lennardJones(double epsilon, double sigma, double cutOff);

        double evaluate(particle &part1, particle &part2) override;

        std::array<vec3<double>, 4>
        forceTorque(particle &part1, particle &part2) override;
    };

}
