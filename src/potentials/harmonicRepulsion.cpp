//
// Created by maojrs on 8/16/18.
//
//#include <cmath>
#include "vec3.hpp"
#include "potentials/harmonicRepulsion.hpp"

/** Class implementation for harmonic repulsion pair potential
 * @param k WHAT IS k?
 * @param range WHAT IS range?
 */
harmonicRepulsion::harmonicRepulsion(double k, double range) : k(k), range(range) {}


double harmonicRepulsion::evaluate(vec3<double> pos1, vec3<double> pos2) {
    vec3<double> d = pos2 - pos1;
    double R = d.norm();
    if (R > range) {
        return 0;
    } else {
        return k/2*std::pow(R-range, 2);
    }
}

// Returns -gradient of potential (force) at position x
vec3<double> harmonicRepulsion::force(vec3<double> pos1, vec3<double> pos2) {
    vec3<double> d = pos2-pos1;
    double R = d.norm();
    if (R > range) {
        return 0.*d;
    } else {
        return k*(R-range) / R * d;
    }
};

