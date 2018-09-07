//
// Created by maojrs on 8/21/18.
//

#pragma once
#include "potentials.hpp"

namespace msmrd {
    /*
     * Declaration of dipole potential (models dipole under constant isotropic electric field)
     * The template <vec3<double>> indicates potential depends on orientation described by a vector (rod-like particles)
     */
    class dipole : public externalPotential<vec3<double>> {
    public:
        double scalefactor = 1.0;
        vec3<double> directionEField = vec3<double>(0, 0, 1);

        /**
        * @param scalefactor adjusts strength of the potential muliplying by a constant
        * @param directionEField represents the direction of the electric field
        */
        dipole(double scalefactor, vec3<double> directionEField);

        dipole(double scalefactor, std::vector<double> &directionEField);


        double evaluate(vec3<double> pos, vec3<double> u) override;

        std::array<vec3<double>, 2> forceTorque(vec3<double> pos, vec3<double> u) override;
    };

}
