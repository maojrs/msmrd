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

double harmonicRepulsion::evaluatePyBind(std::vector<double> pos1, std::vector<double> pos2) {
    vec3<double> x = vec3<double>(pos1);
    vec3<double> y = vec3<double>(pos2);
    return evaluate(x, y);
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

// Needed for PyBinding force, since it can take vectors as input.
std::vector<double> harmonicRepulsion::forcePyBind(std::vector<double> pos1, std::vector<double> pos2) {
    vec3<double> x = vec3<double>(pos1);
    vec3<double> y = vec3<double>(pos2);
    vec3<double> forcex = force(x, y);
    std::vector<double> output;
    output[0] = forcex[0];
    output[1] = forcex[1];
    output[2] = forcex[2];
    return output;
}

