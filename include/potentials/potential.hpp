//
// Created by maojrs on 8/14/18.
//
#pragma once
#include "vec3.hpp"
#include "randomgen.hpp"

/**
 * Abstract base class definition for external potentials
 */
class potential{
public:
    potential() = default;

     // Calculate value of potential and force at position "pos"
    virtual double evaluate(vec3<double> pos) = 0;
    virtual vec3<double> force(vec3<double> pos) = 0;
};


/**
 * Abstract base class definition for pair potentials
 */
class pairPotential{
public:
    pairPotential() = default;

    // Calculate value of potential and force at position "pos"
    virtual double evaluate(vec3<double> pos1, vec3<double> pos2) = 0;
    virtual vec3<double> force(vec3<double> pos1, vec3<double> pos2) = 0;
};


// Defines and implements null external potential
class nullPotential: public potential {
public:
    double evaluate(vec3<double> pos) {return 0;}
    vec3<double> force(vec3<double> pos) {return vec3<double>(0,0,0);}
};


// Defines and implements null pairPotential
class nullPairPotential: public pairPotential {
public:
    double evaluate(vec3<double> pos1, vec3<double> pos2) {return 0;}
    vec3<double> force(vec3<double> pos1, vec3<double> pos2) {return vec3<double>(0,0,0);}
};




