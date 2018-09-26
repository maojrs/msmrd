//
// Created by maojrs on 9/26/18.
//

#pragma once
#include "potentials/patchyParticle.hpp"

namespace msmrd{

    /*
     * Declaration of potential function for patchy protein model. The class is a child from the
     * general potential for patchy particles. Note the potential depends on the
     * orientation of both particles, each with all its three rotational degrees of freedom.
     */
    class patchyProtein : public patchyParticle {
    public:
        patchyProtein(double sigma, double strength, std::vector<vec3<double>> patchesCoordinates);
        patchyProtein(double sigma, double strength, std::vector<std::vector> patchesCoordinates);

        double evaluate(vec3<double> pos1, vec3<double> pos2, quaternion<double> theta1, quaternion<double> theta2) override;

        std::array<vec3<double>, 4>
        forceTorque(vec3<double> pos1, vec3<double> pos2, quaternion<double> theta1, quaternion<double> theta2) override;

    };


}

