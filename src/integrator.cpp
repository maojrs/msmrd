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
        integrateOne(parts[i]);
    }
    clock += dt;
};

/**
 * Functions and constructors of child classes of integrator
 */

// Over-damped Lanegvin dynamics integrator
odLangevin::odLangevin(double dt, long seed, bool rotation) : integrator(dt,seed, rotation) {};

void odLangevin::integrateOne(particle &part) {
    translate(part,dt);
    if (rotation) {
        rotate(part, dt);
    }
}

void odLangevin::integrate(particle &part) {
    integrateOne(part);
    clock += dt;
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

/**
* Over-damped Langevin dynamics with Markovian switch integrator
* constructors defined in headers since they use templates.
*/

 
// Integrates rotation/translation and Markovian switch together
template<>
void odLangevinMarkovSwitch<ctmsm>::integrateOne(particleMS &part) {
    double resdt;
    // propagate CTMSM/MSM when synchronized and update diffusion coefficients
    part.tcount = 0;
    // Runs one timestep dt, if lagtime < dt, propagates MSM as many times as needed
    if (part.lagtime <= dt) {
        // Run loop until integration for one timestep dt is achieved
        while (part.tcount < dt) {
            // Propagates MSM only when diffusion and rotation are in sync
            if (part.propagateTMSM) {
                tmsm.propagate(part, 1);
                part.setD(tmsm.Dlist[part.state]);
                part.setDrot(tmsm.Drotlist[part.state]);
            }
            // Integrates for one lagtime as long as integration is still under dt
            if ( part.tcount + part.lagtime < dt ) {
                translate(part, part.lagtime);
                rotate(part, part.lagtime);
                part.tcount += part.lagtime;
                part.propagateTMSM = true;
            // If current lagtime overtakes dt, integrate up to dt (by resdt) and reset lagtime to remaining portion
            } else {
                resdt = dt - part.tcount;
                translate(part, resdt);
                rotate(part, resdt);
                part.setLagtime(part.lagtime + part.tcount - dt);
                part.tcount += resdt; // this means part.tcount = dt, so will exit while loop.
                // If lag time = 0 MSM must propagate in next step, otherwise it needs to integrate remaining lagtime.
                if (part.lagtime == 0) {
                    part.propagateTMSM = true;
                } else {
                    part.propagateTMSM = false;
                };
            };
        }
        part.tcount = 0;
    } else {
    // Runs one full time step when lagtime > dt and update remaining lagtime.
        translate(part, dt);
        rotate(part, dt);
        part.setLagtime(part.lagtime - dt);
        // If lag time = 0 MSM must propagate in next step, otherwise it needs to integrate remaining lagtime.
        if (part.lagtime == 0) {
            part.propagateTMSM = true;
        } else {
            part.propagateTMSM = false;
        };
    };
};

template<>
void odLangevinMarkovSwitch<ctmsm>::integrate(particleMS &part) {
    integrateOne(part);
    clock += dt;
};



// Integrate list of particlesMS (override parent function)
template<>
void odLangevinMarkovSwitch<ctmsm>::integrateList(std::vector<particleMS> &parts) {
    for (int i=0; i<parts.size(); i++) {
        integrateOne(parts[i]);
    }
    clock += dt;
};



