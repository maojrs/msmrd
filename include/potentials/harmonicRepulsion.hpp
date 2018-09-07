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
    class harmonicRepulsion : public pairPotential<> {
    public:
        double k;
        double range;

        /**
         * @param k repulsion strength
         * @param range interaction radius
         */
        harmonicRepulsion(double k, double range);

        double evaluate(vec3<double> pos1, vec3<double> pos2) override;

        std::array<vec3<double>, 2> forceTorque(vec3<double> pos1, vec3<double> pos2) override;
    };

}