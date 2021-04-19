//
// Created by maojrs on 4/19/21.
//

#include "potentials/lennardJones.hpp"

namespace msmrd {
    /**
     * Class implementation for Lennard Jonnes pair potential
     */
    lennardJones::lennardJones(double epsilon, double sigma) : epsilon(epsilon), sigma(sigma) {}

    double lennardJones::evaluate(particle &part1, particle &part2) {
        std::array<vec3<double>, 2> relPos = relativePositionComplete(part1.position, part2.position);
        //vec3<double> pos1virtual = relPos[0]; // virtual pos1 if periodic boundary; otherwise pos1.
        vec3<double> rvec = relPos[1]; //part2.position - part1.position if non-periodic boundary;
        double r = rvec.norm();
        if (r <= 3*sigma) {
            auto term0 = std::pow((sigma/r),6);
            auto term1 = std::pow(term0,2);
            return 4*epsilon*(term1 - term0);
        } else {
            return 0;
        }
    }

    std::array<vec3<double>, 4> lennardJones::forceTorque(particle &part1, particle &part2) {
    // MISSING IMPLEMENTATION
    }


}