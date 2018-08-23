//
// Created by maojrs on 8/21/18.
//
#include "potentials/dipole.hpp"

/**
 * Class implementation for dipole potential (dipole under constant isotropic electric field)
 * @param scalefactor adjusts strength of the potential muliplying by a constant
 * @param directionEField represents the direction of the electric field
 */
dipole::dipole(double scalefactor, vec3<double> directionEField)
        : scalefactor(scalefactor), directionEField(directionEField) {};

dipole::dipole(double scalefactor, std::vector<double> &directionEField)
        : scalefactor(scalefactor), directionEField(directionEField) {};

double dipole::evaluate(vec3<double> pos, vec3<double> u) {
 return pos*u;
};

std::array<vec3<double>, 2> dipole::forceTorque(vec3<double> pos, vec3<double> u) {
    vec3<double> force = vec3<double>(0., 0., 0.);
    vec3<double> torque;
    torque = u.cross(directionEField);
    return {force, scalefactor*torque};
};