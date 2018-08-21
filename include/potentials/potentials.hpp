//
// Created by maojrs on 8/14/18.
//
#pragma once
#include "vec3.hpp"
#include "randomgen.hpp"

/**
 * Abstract base class declaration for external potentials
 */
class externalPotential{
public:
    externalPotential() = default;

     // Calculate value of potential and force at position "pos"
    virtual double evaluate(vec3<double> pos) = 0;
    double evaluatePyBind(std::vector<double> pos);
    virtual vec3<double> force(vec3<double> pos) = 0;
    std::vector<double> forcePyBind(std::vector<double> pos);
};


/**
 * Abstract base class declaration for external potentials for rod-like particles
 */
class externalRodPotential{
public:
    externalRodPotential() = default;

    // Calculate value of potential and force at position "pos"
    virtual double evaluate(vec3<double> pos, vec3<double> u) = 0;
    double evaluatePyBind(std::vector<double> pos, std::vector<double> u);
    virtual std::array<vec3<double>, 2> forceTorque(vec3<double> pos, vec3<double> u) = 0;
    std::vector<std::vector<double>> forcePyBind(std::vector<double> pos, std::vector<double> u);
};


/**
 * Abstract base class declaration for pair potentials
 */
class pairPotential{
public:
    pairPotential() = default;

    // Calculate value of potential and force at position "pos"
    virtual double evaluate(vec3<double> pos1, vec3<double> pos2) = 0;
    double evaluatePyBind(std::vector<double> pos1, std::vector<double> pos2);
    virtual vec3<double> force(vec3<double> pos1, vec3<double> pos2) = 0;
    std::vector<double> forcePyBind(std::vector<double> pos1, std::vector<double> pos2);
};


/**
 * Abstract base class declaration for pair potentials between rod-like particles,
 * which means their orientation can be determined from one vector
 */
class rodPairPotential{
public:
    rodPairPotential() = default;
    // Calculate value of potential and force at position "pos"
    virtual double evaluate(vec3<double> pos1, vec3<double> pos2, vec3<double> u1, vec3<double> u2) = 0;
    double evaluatePyBind(std::vector<double> pos1, std::vector<double> pos2,std::vector<double> u1, std::vector<double> u2);
    virtual std::array<vec3<double>, 2> forceTorque(vec3<double> pos1, vec3<double> pos2, vec3<double> u1, vec3<double> u2) = 0;
    std::vector<std::vector<double>> forcePyBind(std::vector<double> pos1, std::vector<double> pos2,std::vector<double> u1, std::vector<double> u2);
};


/**
 *  Null potentials for abstract potentials classes declaration and implementation
 */

// Defines and implements null external potential
class nullExternalPotential: public externalPotential {
public:
    double evaluate(vec3<double> pos) {return 0;}
    vec3<double> force(vec3<double> pos) {return vec3<double>(0,0,0);}
};

// Defines and implements null external potential
class nullExternalRodPotential: public externalPotential {
public:
    double evaluate(vec3<double> pos, vec3<double> u) {return 0;}
    vec3<double> force(vec3<double> pos, vec3<double> u) {return vec3<double>(0,0,0);}
};

// Defines and implements null pairPotential
class nullPairPotential: public pairPotential {
public:
    double evaluate(vec3<double> pos1, vec3<double> pos2) {return 0;}
    vec3<double> force(vec3<double> pos1, vec3<double> pos2) {return vec3<double>(0,0,0);}
};

// Defines and implements null pairPotentialTorque
class nullRodPairPotential: public pairPotential {
public:
    double evaluate(vec3<double> pos1, vec3<double> pos2, vec3<double> u1, vec3<double> u2) {return 0;}
    std::array<vec3<double>, 2> forceTorque(vec3<double> pos1, vec3<double> pos2, vec3<double> u1, vec3<double> u2) {
        return std::array<vec3<double>, 2>{vec3<double>(0,0,0), vec3<double>(0,0,0)};}
};







