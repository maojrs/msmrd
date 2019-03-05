//
// Created by maojrs on 8/16/18.
//
//#include <cmath>
#include "vec3.hpp"
#include "potentials/harmonicRepulsion.hpp"

namespace msmrd {
    /** Class implementation for harmonic repulsion pair potential
     * @param k repulsion strength
     * @param range interaction radius
     */
    harmonicRepulsion::harmonicRepulsion(double k, double range) : k(k), range(range) {}

    // Evaluate potential value for two given particles' positions
    double harmonicRepulsion::evaluate(const particle &part1, const particle &part2) {
        vec3<double> d = relativePosition(part2.position, part1.position); //part1.position - part2.position;
        double R = d.norm();
        if (R > range) {
            return 0;
        } else {
            return k / 2 * std::pow(R - range, 2);
        }
    }

    // Returns -gradient of potential (force)  and zero torque at position x
    std::array<vec3<double>, 4> harmonicRepulsion::forceTorque(const particle &part1, const particle &part2) {
        vec3<double> force = vec3<double>(0, 0, 0);
        vec3<double> torque = vec3<double>(0, 0, 0);
        vec3<double> d = relativePosition(part2.position, part1.position); //part1.position - part2.position;
        double R = d.norm();
        if (R > range) {
            force = 0. * d;
        } else {
            force = k * (R - range) / R * d;
        }
        return {force, torque, -1.0*force, -1.0*torque};
    }

}

