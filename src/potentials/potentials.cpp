//
// Created by maojrs on 8/21/18.
//
#include <tuple>
#include "potentials/potentials.hpp"
#include "tools.hpp"
#include "vec3.hpp"

namespace msmrd {
    /**
     *  Implementation of non-abstract functions of externalPotential abstract class (default constructor in header)
     *  These functions are needed for PyBinding forceTorque functions, since pyBind only vector output.
     */


    /* Force Torque for binding functions */
    std::vector<std::vector<double>> externalPotential::forceTorquePyBind(const particle &part) {
        std::array<vec3<double>, 2> forceTorquex = forceTorque(part);
        return msmrdtools::array2Dtovec2D(forceTorquex);
    }

    std::vector<std::vector<double>> pairPotential::forceTorquePyBind(const particle &part1, const particle &part2) {
        std::array<vec3<double>, 4> forceTorquex = forceTorque(part1, part2);
        return msmrdtools::array2Dtovec2D(forceTorquex);
    }

}