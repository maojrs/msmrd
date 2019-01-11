//
// Created by maojrs on 8/21/18.
//

#pragma once

#include "particle.hpp"
#include "potentials.hpp"

namespace msmrd {
    /*
     * Declaration of dipole potential (models dipole under constant isotropic electric field)
     * The template <vec3<double>> indicates potential depends on orientation described by a vector (rod-like particles)
     */
    class dipole : public externalPotential {
    public:
        double scalefactor = 1.0;
        vec3<double> directionEField = vec3<double>(0, 0, 1);

        /**
        * @param scalefactor adjusts strength of the potential muliplying by a constant
        * @param directionEField represents the direction of the electric field
        */
        dipole(double scalefactor, vec3<double> directionEField);

        dipole(double scalefactor, std::vector<double> &directionEField);


        double evaluate(const particle &part) override;

        std::array<vec3<double>, 2> forceTorque(const particle &part) override;
    };

}
