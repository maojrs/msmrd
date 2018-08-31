//
// Created by maojrs on 8/16/18.
//

#include "integrators/odLangevinMarkovSwitch.hpp"
#include "particle.hpp"
#include "msm.hpp"

/**
* Implementation of over-damped Langevin dynamics with Markovian switch integrator class
* (constructors defined in headers since templates were used).
*/


// Integrates rotation/translation and Markovian switch of one particle (visible only inside the class)
template<>
void odLangevinMarkovSwitch<ctmsm>::integrateOne(particleMS &part) {
    // Calculate forces and torque in this time step
    vec3<double> force;
    vec3<double> torque;
    std::array<vec3<double>, 2> forctorq;
    // Do diffusion propagation taking MSM/CTMSM into account
    double resdt;
    // propagate CTMSM/MSM when synchronized and update diffusion coefficients
    part.tcount = 0;
    // Runs one timestep dt, if lagtime < dt, propagates MSM as many times as needed
    if (part.lagtime <= dt) {
        // Run loop until integration for one timestep dt is achieved
        while (part.tcount < dt) {
            // Propagates MSM only when diffusion and rotation are in sync
            if (part.propagateTMSM) {
                tmsm.propagate(part, 1); // Diffusion coefficients don't need to be updated until dt reaches lagtime.
            }
            // Integrates for one lagtime as long as integration is still under dt
            if ( part.tcount + part.lagtime < dt ) {
                forctorq = getExternalForceTorque(part);
                force = forctorq[0];
                torque = forctorq[1];
                translate(part, force, part.lagtime);
                rotate(part, torque, part.lagtime);
                part.tcount += part.lagtime;
                part.setLagtime(0);
                // Ready to propagate MSM, update state and diffusion coefficients
                part.propagateTMSM = true;
                part.setState(part.nextState);
                part.setD(tmsm.Dlist[part.state]);
                part.setDrot(tmsm.Drotlist[part.state]);
                // If current lagtime overtakes dt, integrate up to dt (by resdt) and reset lagtime to remaining portion
            } else {
                forctorq = getExternalForceTorque(part);
                force = forctorq[0];
                torque = forctorq[1];
                resdt = dt - part.tcount;
                translate(part, force, resdt);
                rotate(part, torque, resdt);
                part.setLagtime(part.lagtime + part.tcount - dt);
                part.tcount += resdt; // this means part.tcount = dt, so will exit while loop.
                // If lag time = 0 MSM must propagate in next step, otherwise it needs to integrate remaining lagtime.
                if (part.lagtime == 0) {
                    // Ready to propagate MSM, update state and diffusion coefficients
                    part.propagateTMSM = true;
                    part.setState(part.nextState);
                    part.setD(tmsm.Dlist[part.state]);
                    part.setDrot(tmsm.Drotlist[part.state]);
                } else {
                    part.propagateTMSM = false;
                };
            };
        }
        part.tcount = 0;
    } else {
        // Runs one full time step when lagtime > dt and update remaining lagtime.
        forctorq = getExternalForceTorque(part);
        force = forctorq[0];
        torque = forctorq[1];
        translate(part, force, dt);
        rotate(part, torque, dt);
        part.setLagtime(part.lagtime - dt);
        // If lag time = 0 MSM must propagate in next step, otherwise it needs to integrate remaining lagtime.
        if (part.lagtime == 0) {
            // Ready to propagate MSM, update state and diffusion coefficients
            part.propagateTMSM = true;
            part.setState(part.nextState);
            part.setD(tmsm.Dlist[part.state]);
            part.setDrot(tmsm.Drotlist[part.state]);
        } else {
            part.propagateTMSM = false;
        };
    };
};


// Same as integrateOne, but updates time and publicly visible
template<>
void odLangevinMarkovSwitch<ctmsm>::integrate(particleMS &part) {
    integrateOne(part);
    clock += dt;
};


// Integrates list of particlesMS (needs to override parent function because it is template based and uses particleMS)
template<>
void odLangevinMarkovSwitch<ctmsm>::integrateList(std::vector<particleMS> &parts) {
    for (int i=0; i<parts.size(); i++) {
        integrateOne(parts[i]);
    }
    clock += dt;
};

