//
// Created by maojrs on 8/16/18.
//
//#include <cmath>
#include "vec3.hpp"
#include "potentials/harmonicRepulsion.hpp"

/** Class implementation for harmonic repulsion pair potential
 * @param k repulsion strength
 * @param range interaction radius
 */
harmonicRepulsion::harmonicRepulsion(double k, double range) : k(k), range(range) {}

// Evaluate potential value for two given particles' positions
double harmonicRepulsion::evaluate(vec3<double> pos1, vec3<double> pos2) {
    vec3<double> d = pos2 - pos1;
    double R = d.norm();
    if (R > range) {
        return 0;
    } else {
        return k/2*std::pow(R-range, 2);
    }
}

// Returns -gradient of potential (force)  and zero torque at position x
std::array<vec3<double>, 2> harmonicRepulsion::forceTorque(vec3<double> pos1, vec3<double> pos2) {
    vec3<double> force = vec3<double>(0, 0, 0);
    vec3<double> torque = vec3<double>(0, 0, 0);
    vec3<double> d = pos2-pos1;
    double R = d.norm();
    if (R > range) {
        force = 0.*d;
    } else {
        force = k*(R-range) / R * d;
    }
    return {force, torque};
};

