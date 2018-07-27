//
// Created by dibakma on 27.06.18.
//
#include <array>
#include <utility>
#include "integrator.hpp"
#include "particle.hpp"
#include "msm.hpp"


/**
 * Functions of parent class integrator
 */
 // Integrate list of particles instead of single one (need to be overriden for interacting particles)
void integrator::integrateList(std::vector<particle> &parts) {
    vec3<double> dr;
    double coeff;
    for (int i=0; i<parts.size(); i++) {
        integrate(parts[i]);
    }
};


/**
 * Functions and constructors of child classes of integrator
 */

// Over-damped Lanegvin dynamics integrator
odLangevin::odLangevin(double dt, long seed) : integrator(dt,seed) {};

void odLangevin::integrate(particle &part) {
    vec3<double> dr;
    double coeff;
    coeff = std::sqrt(2*dt*part.D);
    dr[0] = coeff*randg.normal(0,1);
    dr[1] = coeff*randg.normal(0,1);
    dr[2] = coeff*randg.normal(0,1);
    part.setPosition(part.position + dr);
}

void odLangevin::test(std::vector<int> &intlist){
    intlist[0] = 666;
    //return randg.normal(0,1);
}


// Over-damped Langevin dynamics with Markovian switch integrator
odLangevinMarkovSwitch::odLangevinMarkovSwitch(std::vector<ctmsm> &msmlist, double dt, long seed)
        : msmlist(msmlist), integrator(dt,seed) {};

void odLangevinMarkovSwitch::integrate(particle &part) {
    vec3<double> dr;
    double coeff;
    coeff = std::sqrt(2*dt*part.D);
    dr[0] = coeff*randg.normal(0,1);
    dr[1] = coeff*randg.normal(0,1);
    dr[2] = coeff*randg.normal(0,1);
    part.setPosition(part.position + dr);
};



