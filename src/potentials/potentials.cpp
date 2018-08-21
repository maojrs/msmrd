//
// Created by maojrs on 8/21/18.
//
#include "vec3.hpp"
#include "potentials/potentials.hpp"

/**
 *  Implementation of non-abstract functions of potentials abstract class (default constructor in header)
 */

// Needed for PyBinding externalPotential.evaluate function, since it can only take vectors as input.
double externalPotential::evaluatePyBind(std::vector<double> pos) {
    vec3<double> x = vec3<double>(pos);
    return evaluate(x);
}

// Needed for PyBinding externalPotential.force function , since it can only take vectors as input.
std::vector<double> externalPotential::forcePyBind(std::vector<double> pos) {
    vec3<double> x = vec3<double>(pos);
    vec3<double> forcex = force(x);
    std::vector<double> output;
    output[0] = forcex[0];
    output[1] = forcex[1];
    output[2] = forcex[2];
    return output;
}

// Needed for PyBinding externalRodPotential.evaluate function, since it can only take vectors as input.
double externalRodPotential::evaluatePyBind(std::vector<double> pos1, std::vector<double> u) {
    vec3<double> x = vec3<double>(pos1);
    vec3<double> theta = vec3<double>(u);
    return evaluate(x,theta);
}

// Needed for PyBinding externalRodPotential.force function , since it can only take vectors as input.
std::vector<std::vector<double>> externalRodPotential::forcePyBind(std::vector<double> pos1, std::vector<double> u) {
    vec3<double> x = vec3<double>(pos1);
    vec3<double> theta = vec3<double>(u);
    std::array<vec3<double>, 2> forceTorquex = forceTorque(x, theta);
    std::vector<std::vector<double>> output;
    output.resize(2);
    output[0].resize(3);
    output[1].resize(3);
    for (int i = 0; i<2; i++) {
        output[i][0] = forceTorquex[i][0];
        output[i][1] = forceTorquex[i][1];
        output[i][2] = forceTorquex[i][2];
    }
    return output;
}

// Needed for PyBinding pairPotential.evaluate function, since it can only take vectors as input.
double pairPotential::evaluatePyBind(std::vector<double> pos1, std::vector<double> pos2) {
    vec3<double> x1 = vec3<double>(pos1);
    vec3<double> x2 = vec3<double>(pos2);
    return evaluate(x1,x2);
}

// Needed for PyBinding pairPotential.force function , since it can only take vectors as input.
std::vector<double> pairPotential::forcePyBind(std::vector<double> pos1, std::vector<double> pos2) {
    vec3<double> x1 = vec3<double>(pos1);
    vec3<double> x2 = vec3<double>(pos2);
    vec3<double> forcex = force(x1,x2);
    std::vector<double> output;
    output[0] = forcex[0];
    output[1] = forcex[1];
    output[2] = forcex[2];
    return output;
}

// Needed for PyBinding rodPairPotential.evaluate function, since it can only take vectors as input.
double rodPairPotential::evaluatePyBind(std::vector<double> pos1, std::vector<double> pos2, std::vector<double> u1, std::vector<double> u2) {
    vec3<double> x1 = vec3<double>(pos1);
    vec3<double> x2 = vec3<double>(pos2);
    vec3<double> theta1 = vec3<double>(u1);
    vec3<double> theta2 = vec3<double>(u2);
    return evaluate(x1,x2,theta1,theta2);
}

// Needed for PyBinding rodPairPotential.force function , since it can only take vectors as input.
std::vector<std::vector<double>> rodPairPotential::forcePyBind(std::vector<double> pos1, std::vector<double> pos2, std::vector<double> u1, std::vector<double> u2) {
    vec3<double> x1 = vec3<double>(pos1);
    vec3<double> x2 = vec3<double>(pos2);
    vec3<double> theta1 = vec3<double>(u1);
    vec3<double> theta2 = vec3<double>(u2);
    std::array<vec3<double>, 2> forceTorquex = forceTorque(x1, x2, theta1, theta2);
    std::vector<std::vector<double>> output;
    output.resize(2);
    output[0].resize(3);
    output[1].resize(3);
    for (int i = 0; i<2; i++) {
        output[i][0] = forceTorquex[i][0];
        output[i][1] = forceTorquex[i][1];
        output[i][2] = forceTorquex[i][2];
    }
    return output;
}