//
// Created by maojrs on 8/16/18.
//

#pragma once
#include "vec3.hpp"
#include "potentials/potentials.hpp"

namespace msmrd {
    /*
     * Harmonic repulsion pair potential declaration
     * The empty template <> indicates orientation is not taken into account by this potential.
     */
    class harmonicRepulsion : public pairPotential {
    public:
        double k;
        double range;

        /**
         * @param k repulsion strength
         * @param range interaction radius
         */
        harmonicRepulsion(double k, double range);

        double evaluate(particle &part1, particle &part2) override;

        std::array<vec3<double>, 4> forceTorque(particle &part1, particle &part2) override;
    };

}