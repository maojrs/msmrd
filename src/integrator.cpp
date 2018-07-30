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
    translate(part,dt);
    rotate(part,dt);
}

void odLangevin::translate(particle &part, double dt0){
    vec3<double> dr;
    dr = std::sqrt(2*dt0*part.D)*randg.normal3D(0,1);
    part.setPosition(part.position + dr);
}

void odLangevin::rotate(particle &part, double dt0){
    vec3<double> dphi;
    quaternion<double> dquat;
    dphi = std::sqrt(2*dt0*part.Drot)*randg.normal3D(0,1);
    dquat = angle2quaternion(dphi);
    dquat = dquat * part.orientation;
    part.setOrientation(dquat);
}

double odLangevin::test(particle &part){
    vec3<double> dphi;
    quaternion<double> dquat;
    dphi = std::sqrt(2*part.Drot)*randg.normal3D(0,1);
    dquat = angle2quaternion(dphi);
    //dquat = dquat * part.orientation;
    part.setOrientation(dphi);
    return dphi[1];
}


// Over-damped Langevin dynamics with Markovian switch integrator
// Constructors to define template specializations for ctmsm and msm.
template<typename TMSM>
odLangevinMarkovSwitch<TMSM>::odLangevinMarkovSwitch(msm &msm0, double dt, long seed)
        : tmsm(msm0), odLangevin(dt,seed) {
    msmtype = "discrete-time";
};
template<typename TMSM>
odLangevinMarkovSwitch<TMSM>::odLangevinMarkovSwitch(ctmsm &ctmsm0, double dt, long seed)
        : tmsm(ctmsm0), odLangevin(dt,seed) {
    msmtype = "continuous-time";
};

template<>
void odLangevinMarkovSwitch<ctmsm>::integrate(particle &part) {
    translate(part, dt);
    rotate(part,dt);
}



