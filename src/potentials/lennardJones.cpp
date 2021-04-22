//
// Created by maojrs on 4/19/21.
//

#include "potentials/lennardJones.hpp"

namespace msmrd {
    /**
     * Class implementation for Lennard Jonnes pair potential
     */
    lennardJones::lennardJones(double epsilon, double sigma) : epsilon(epsilon), sigma(sigma) {
        cutOff = 3 * sigma;
    }

    lennardJones::lennardJones(double epsilon, double sigma, double cutOff) : epsilon(epsilon), sigma(sigma),
        cutOff(cutOff) {}

    double lennardJones::evaluate(particle &part1, particle &part2) {
        std::array<vec3<double>, 2> relPos = relativePositionComplete(part1.position, part2.position);
        //vec3<double> pos1virtual = relPos[0]; // virtual pos1 if periodic boundary; otherwise pos1.
        vec3<double> rvec = relPos[1]; //part2.position - part1.position if non-periodic boundary;
        double r = rvec.norm();
        if (r <= cutOff) {
            auto term0 = std::pow((sigma/r),6);
            auto term1 = std::pow(term0,2);
            return 4 * epsilon * (term1 - term0) + baseEnergy;
        } else {
            return 0;
        }
    }

    std::array<vec3<double>, 4> lennardJones::forceTorque(particle &part1, particle &part2) {
        std::array<vec3<double>, 2> relPos = relativePositionComplete(part1.position, part2.position);
        //vec3<double> pos1virtual = relPos[0]; // virtual pos1 if periodic boundary; otherwise pos1.
        vec3<double> rvec = relPos[1]; //part2.position - part1.position if non-periodic boundary;
        double r = rvec.norm();
        vec3<double> force = vec3<double>();
        vec3<double> torque = vec3<double>();
        if (r <= cutOff) {
            double term0 = std::pow((sigma / r), 6);
            double dVdr = (24 * epsilon / r) * (term0 - 2 * std::pow(term0, 2));
            double dVdx = dVdr * rvec[0] / r;
            double dVdy = dVdr * rvec[1] / r;
            double dVdz = dVdr * rvec[2] / r;

            force = vec3<double>(dVdx, dVdy, dVdz);
            torque = 0 * force;
            return {force, torque, -1 * force, torque};
        }
        return {force, torque, -1 * force, torque};
    }


    /**
     * Class implementation for WCA pair potential
     */
    WCA::WCA(double epsilon, double sigma) : lennardJones(epsilon,sigma) {
        cutOff = std::pow(2,1.0/6.0) * sigma;
        baseEnergy = 1 * epsilon;
    }


}