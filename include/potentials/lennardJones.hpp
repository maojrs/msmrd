//
// Created by maojrs on 4/19/21.
//

#pragma once
#include "potentials.hpp"

namespace msmrd {
    /*
     * Declaration of the Lennard Jones pair potential. We also include the WCA potential here, which is simply a
     * Lennard-Jones pontetical truncated at its minima, so it only has the repulsive part.
     */
    class lennardJones : public pairPotential {
    public:
        double epsilon;
        double sigma;
        double cutOff;
        double baseEnergy = 0.0;

        /**
         * @param epsilon depth of the potential well
         * @param sigma distance at which the potential is zero (a.k.a 'size of the particle')
         * @param cutOff value after which potentil is zero, defaults to 3 * sigma if not specified
         * @param baseEnergy translation of potential energy in y-axis. It defaults to zero, but it can
         * be modified by child classes, like the WCA.
         */
        lennardJones(double epsilon, double sigma);

        lennardJones(double epsilon, double sigma, double cutOff);

        double evaluate(particle &part1, particle &part2) override;

        std::array<vec3<double>, 4>
        forceTorque(particle &part1, particle &part2) override;
    };

    class WCA : public lennardJones {
    public:
        WCA(double epsilon, double sigma);
    };

}
