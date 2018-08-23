//
// Created by dibakma on 27.06.18.
//
#include <array>
#include <utility>
#include "integrators/integrator.hpp"
#include "particle.hpp"


/**
 * Implementation of integrator abstract class inherited by all child classes
 */
integrator::integrator(double dt, long seed, bool rotation)
        : dt(dt), seed(seed), rotation(rotation) {
    // Define zero potentials
    nullExternalPotential nullpot;
    nullExternalRodPotential nullrodpot;
    nullPairPotential nullpairpot;
    nullRodPairPotential nullpairrodpot;
    // set all potential types to zero as default
    externalPot = &nullpot;
    externalRodPot = &nullrodpot;
    pairPot = &nullpairpot;
    pairRodPot = &nullpairrodpot;
    randg.setSeed(seed);
    clock = 0;
};

// Calculates force and torque from an external potential for point, rod-like and rigidsolid particles
std::array<vec3<double>, 2> integrator::getExternalForceTorque(particle &part) {
    if (part.bodytype == "point") {
        return externalPot->forceTorque(part.position);
    } else if (part.bodytype == "rod") {
        return externalRodPot->forceTorque(part.position, part.orientvector);
    } else if (part.bodytype == "rigidsolid") {
        // to be implemented (forceTorque functions in potentials from quaternions)
        // return externalQuatPot->forceTorque(part.position, part.orientation);
    } else {
        throw std::runtime_error("Unknown particle bodytype. it should be either point, rod or rigidsolid.");
    };
};

// Calculates force and torque due to pair interactions for point, rod-like and rigidsolid particles
std::array<vec3<double>, 2> integrator::getPairsForceTorque(particle &part, std::vector<particle> &parts) {
    std::array<vec3<double>, 2> forctorq;
    vec3<double> force = vec3<double>(0.,0.,0.);
    vec3<double> torque = vec3<double>(0.,0.,0.);
    if (part.bodytype == "point") {
        for (int i=0; i<parts.size(); i++) {
            forctorq = pairPot->forceTorque(part.position, parts[i].position);
            force += forctorq[0];
            torque += forctorq[1];
        }
        return {force, torque};
    } else if (part.bodytype == "rod") {
        for (int i=0; i<parts.size(); i++) {
            forctorq = pairRodPot->forceTorque(part.position, part.orientvector, parts[i].position, parts[i].orientvector);
            force += forctorq[0];
            torque += forctorq[1];
        }
        return {force, torque};
    } else if (part.bodytype == "rigidsolid") {
        // to be implemented (forceTorque functions in potentials from quaternions)
//        for (int i=0; i<parts.size(); i++) {
//            forctorq = quatPairPot->forceTorque(part.position, part.orientvector, parts[i].position, parts[i].orientvector);
//            force += forctorq[0];
//            torque += forctorq[1];
//        }
//        return {force, torque};
    } else {
        throw std::runtime_error("Unknown particle bodytype. it should be either point, rod or rigidsolid.");
    };
};



// Incorporates custom external potential functions into integrator
void integrator::setExternalPotential(externalPotential<> *pot) {  externalPot = pot; }

// Incorporates custom external potential functions into integrator
void integrator::setExternalRodPotential(externalPotential<vec3<double>> *pot) {  externalRodPot = pot; }

// Incorporates custom pair potential functions into integrator
void integrator::setPairPotential(pairPotential<> *pot) { pairPot = pot; }

// Incorporates custom pair potential function for rod-like particles into integrator
void integrator::setRodPairPotential(pairPotential<vec3<double>,vec3<double>> *pot) { pairRodPot = pot; }


//// Evaluates external potential from integrator at a given position
//double integrator::evalExternalPotential(std::vector<double> pos) {
//    return externalPot->evaluatePyBind(pos);
//}
//
//// Evaluates pair potential from integrator at a given position
//double integrator::evalPairPotential(std::vector<double> pos1, std::vector<double> pos2) {
//    return pairPot->evaluatePyBind(pos1, pos2);
//}
//
//// Evaluates pair potential of rod-like particles from integrator at a given position
//double integrator::evalRodPairPotential(std::vector<double> pos1, std::vector<double> pos2,
//                                        std::vector<double> u1, std::vector<double> u2) {
//    return rodPairPot->evaluatePyBind(pos1, pos2, u1, u2);
//}
//
//// Evaluates force due to external potential from integrator at a given position
//std::vector<double> integrator::evalExternalForce(std::vector<double> pos) {
//    return externalPot->forcePyBind(pos);
//}
//
//// Evaluates force due to pair potential from integrator at a given position
//std::vector<double> integrator::evalPairForce(std::vector<double> pos1, std::vector<double> pos2) {
//    return pairPot->forcePyBind(pos1, pos2);
//}
//
//// Evaluates force due to pair potential from integrator at a given position
//std::vector<std::vector<double>> integrator::evalRodPairForce(std::vector<double> pos1, std::vector<double> pos2,
//                                                              std::vector<double> u1, std::vector<double> u2) {
//    return rodPairPot->forcePyBind(pos1, pos2, u1, u2);
//}

 // Integrate list of particles instead of single one (need to override in case of interacting or MS particles)
void integrator::integrateList(std::vector<particle> &parts) {
    for (int i=0; i<parts.size(); i++) {
        integrateOne(parts[i]);
    }
    clock += dt;
}



