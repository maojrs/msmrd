//
// Created by maojrs on 8/21/18.
//
#include "potentials/dipole.hpp"

namespace msmrd {
    /**
     * Class implementation for dipole potential (dipole under constant isotropic electric field)
     * @param scalefactor adjusts strength of the potential muliplying by a constant
     * @param directionEField represents the direction of the electric field
     */
    dipole::dipole(double scalefactor, vec3<double> directionEField)
            : scalefactor(scalefactor), directionEField(directionEField) {};

    dipole::dipole(double scalefactor, std::vector<double> &directionEField)
            : scalefactor(scalefactor), directionEField(directionEField) {};

    double dipole::evaluate(particle &part) {
        return part.position * part.orientvector;
    };

// Calculates force and torque due to potential, output force is zero
    std::array<vec3<double>, 2> dipole::forceTorque(particle &part) {
        vec3<double> force = vec3<double>(0., 0., 0.);
        vec3<double> torque;
        torque = part.orientvector.cross(directionEField);
        return {force, scalefactor * torque};
    }

}