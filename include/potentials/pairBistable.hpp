//
// Created by maojrs on 2/8/22.
//


#pragma once
#include "vec3.hpp"
#include "potentials/potentials.hpp"

namespace msmrd {
    /*
     * Bistable pair potential for dimers with two possible configurations
     * The potential is of the form: scaleFactor ( 1 - ((x-(x0+rad))/rad)^2 )^2
     */
    class pairBistable : public pairPotential {
    public:
        double x0;
        double rad;
        double scalefactor;
        std::vector<int> particleTypes;

        /**
         * @param x0 location of first minima
         * @param rad half the distance between minimas (distance between minimas = 2 * rad)
         * @param scalefactor factor that scales the whole potential
         * @param particleTypes if assigned values, the external potential only acts on the specified
         * particle types.
         */

        pairBistable(double x0, double rad, double scalefactor);

        pairBistable(double x0, double rad, std::vector<int> particleTypes, double scalefactor);

        double evaluate(particle &part1, particle &part2) override;

        std::array<vec3<double>, 4> forceTorque(particle &part1, particle &part2) override;
    };

    /*
     * Bistable pair potential for dimers with two possible configurations with a bias towards the close configuration
     * The potential is of the form: scaleFactor ( 1 - ((x-(x0+rad))/rad)^2 )^2 + a*log(x)
     */
    class pairBistableBias : public pairPotential {
    public:
        double x0;
        double rad;
        double scalefactor;
        double a = 2.0;
        std::vector<int> particleTypes;

        /**
         * @param x0 approx location of first minima
         * @param rad approx half the distance between minimas (distance between minimas = 2 * rad)
         * @param scalefactor factor that scales the whole potential
         * @param particleTypes if assigned values, the external potential only acts on the specified
         * particle types.
         */

        pairBistableBias(double x0, double rad, double scalefactor);

        pairBistableBias(double x0, double rad, std::vector<int> particleTypes, double scalefactor);

        double evaluate(particle &part1, particle &part2) override;

        std::array<vec3<double>, 4> forceTorque(particle &part1, particle &part2) override;
    };

}
