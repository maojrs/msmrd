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
    nullExternalPotential nullpot;
    nullPairPotential nullpairpot;
    // set all potential types to zero as default
    externalPot = &nullpot;
    pairPot = &nullpairpot;
    randg.setSeed(seed);
    clock = 0;
};

// Incorporates custom external potential functions into integrator
void integrator::setExternalPotential(externalPotential *pot) {  externalPot = pot; }

// Incorporates custom pair potential functions into integrator
void integrator::setPairPotential(pairPotential *pot) { pairPot = pot; }

// Incorporates custom pair potential function for rod-like particles into integrator
void integrator::setRodPairPotential(rodPairPotential *pot) { rodPairPot = pot; }


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



