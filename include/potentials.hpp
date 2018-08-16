//
// Created by maojrs on 8/14/18.
//
#pragma once
#include "vec3.hpp"
#include "randomgen.hpp"

// Abstract base class for potentials
class potentials{
public:
    potentials() = default;

     // Calculate value of potential and force at position "pos"
    virtual double evaluate(vec3<double> pos) = 0;
    virtual vec3<double> force(vec3<double> pos) = 0;
};

// Abstract base class for pair potentials
class pairPotentials{
public:
    pairPotentials() = default;

    // Calculate value of potential and force at position "pos"
    virtual double evaluate(vec3<double> pos1, vec3<double> pos2) = 0;
    virtual double evaluatePyBind(std::vector<double> pos1, std::vector<double> pos2) = 0;
    virtual vec3<double> force(vec3<double> pos1, vec3<double> pos2) = 0;
    virtual std::vector<double> forcePyBind(std::vector<double> pos1, std::vector<double> pos2) = 0;
};


// Null potential
class nullpotential: public potentials {
public:
    double evaluate(vec3<double> pos) {return 0;}
    vec3<double> force(vec3<double> pos) {return vec3<double>(0,0,0);}
};

// 3D potential composed of nminima Gaussians placed randomly inside sphere of radius maxrad
class gaussians3D: public potentials {
private:
    randomgen randg;
public:
    int nminima;
    double maxrad;
    double scalefactor;
    long seed;
    std::vector<vec3<double>> minimas;
    std::vector<vec3<double>> stddevs;
    gaussians3D(int nminima, double maxrad, double scalefactor, long seed);

    double evaluate(vec3<double> pos) override;
    double evaluatePyBind(std::vector<double> pos);
    vec3<double> force(vec3<double> pos) override;
    std::vector<double> forcePyBind(std::vector<double> pos);
};

// harmonic repulsion between particles
class harmonicRepulsion: public pairPotentials{
    public:
    double k;
    double range;
    harmonicRepulsion(double k, double range) : k(k), range(range) {};
    double evaluate(vec3<double> pos1, vec3<double> pos2);
    double evaluatePyBind(std::vector<double> pos1, std::vector<double> pos2);
    vec3<double> force(vec3<double> pos1, vec3<double> pos2);
    std::vector<double> forcePyBind(std::vector<double> pos1, std::vector<double> pos2);
};