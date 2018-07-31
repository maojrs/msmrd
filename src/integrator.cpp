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
 // Integrate list of particles instead of single one (need to be overriden for interacting or MS particles)
void integrator::integrateList(std::vector<particle> &parts) {
    for (int i=0; i<parts.size(); i++) {
        integrate(parts[i]);
    }
};

/**
 * Functions and constructors of child classes of integrator
 */

// Over-damped Lanegvin dynamics integrator
odLangevin::odLangevin(double dt, long seed, bool rotation) : integrator(dt,seed, rotation) {};

void odLangevin::integrate(particle &part) {
    translate(part,dt);
    if (rotation) {
        rotate(part, dt);
    }
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
    part.setOrientation(dquat * part.orientation);
}

// Over-damped Langevin dynamics with Markovian switch integrator
// constructors defined in headers since they use templates.

// Integrates rotation/translation and Markovian switch together
template<>
void odLangevinMarkovSwitch<ctmsm>::integrate(particleMS &part) {
    double resdt;
    // propagate CTMSM/MSM when synchronized and update diffusion coefficients
    if (part.propagateTMSM || part.lagtime == 0){
        tmsm.propagate(part, 1);
        part.tcount = 0;
        part.setD(tmsm.Dlist[part.state]);
        part.setDrot(tmsm.Drotlist[part.state]);
    };
    // integrate full timestep only if current time step is smaller than lagtime
    if (part.tcount + dt < part.lagtime ) {
        translate(part, dt);
        rotate(part,dt);
        part.tcount += dt;
        part.propagateTMSM = false;
        clock += dt;
    }
    // integrate residual time step only if next full time step is beyond lagtime
    // this synchronizes the integration and MSM propagation
    else {
        resdt = part.lagtime - part.tcount;
        translate(part, resdt);
        rotate(part, resdt);
        part.propagateTMSM = true;
        clock += resdt;
    };
}

// Integrate list of particlesMS (override parent function)
template<>
void odLangevinMarkovSwitch<ctmsm>::integrateList(std::vector<particleMS> &parts) {
    for (int i=0; i<parts.size(); i++) {
        integrate(parts[i]);
    }
};



