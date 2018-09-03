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

// Integrates diffusion and rotation of one particle, called by integrateOneMS (visible only inside the class)
template<>
void odLangevinMarkovSwitch<ctmsm>::integrateOne(particleMS &part, double timestep){
    vec3<double> force;
    vec3<double> torque;
    std::array<vec3<double>, 2> forctorq;
    forctorq = getExternalForceTorque(part);
    force = forctorq[0];
    torque = forctorq[1];
    translate(part, force, timestep);
    if (rotation) {
        rotate(part, torque, timestep);
    }
};


/* Integrates diffusion and rotation of one particle wit pair interactions, called by the
 *  corresponding integrateOneMS (visible only inside the class) */
template<>
void odLangevinMarkovSwitch<ctmsm>::integrateOne(int partIndex, std::vector<particleMS> &parts, double timestep){
    vec3<double> force;
    vec3<double> torque;
    std::array<vec3<double>, 2> forctorq;
    forctorq = getExternalForceTorque(parts[partIndex]);
    force = forctorq[0];
    torque = forctorq[1];
    translate(parts[partIndex], force, timestep);
    if (rotation) {
        rotate(parts[partIndex], torque, timestep);
    }
};


/* Integrates rotation/translation and Markovian switch of one particle, without
 * pair interactiones (visible only inside the class) */
template<>
void odLangevinMarkovSwitch<ctmsm>::integrateOneMS(particleMS &part, double timestep) {
    // Do diffusion/rotation propagation taking MSM/CTMSM into account
    double resdt;
    // propagate CTMSM when synchronized and update diffusion coefficients
    part.tcount = 0;
    // Runs one timestep dt, if lagtime < dt, propagates CTMSM as many times as needed
    if (part.lagtime <= timestep) {
        // Run loop until integration for one timestep dt is achieved
        while (part.tcount < timestep) {
            // Propagates MSM only when diffusion and rotation are in sync
            if (part.propagateTMSM) {
                tmsm.propagateNoUpdate(part, 1); // Diffusion coefficients don't need to be updated until dt reaches lagtime.
            }
            // Integrates for one lagtime as long as integration is still under dt
            if ( part.tcount + part.lagtime < timestep ) {
                integrateOne(part, part.lagtime);
                part.tcount += part.lagtime;
                part.setLagtime(0);
                // Ready to propagate MSM, update state and diffusion coefficients
                part.propagateTMSM = true;
                part.setState(part.nextState);
                part.setD(tmsm.Dlist[part.state]);
                part.setDrot(tmsm.Drotlist[part.state]);
                // If current lagtime overtakes dt, integrate up to dt (by resdt) and reset lagtime to remaining portion
            } else {
                resdt = timestep - part.tcount;
                integrateOne(part, resdt);
                part.setLagtime(part.lagtime + part.tcount - timestep);
                part.tcount += resdt; // this means part.tcount = dt, so will exit while loop.
                // If lag time = 0 CTMSM must propagate in next step, otherwise it needs to integrate remaining lagtime.
                if (part.lagtime == 0) {
                    // Ready to propagate CTMSM, update state and diffusion coefficients
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
        integrateOne(part, timestep);
        part.setLagtime(part.lagtime - timestep);
        // If lag time = 0 MSM must propagate in next step, otherwise it needs to integrate remaining lagtime.
        if (part.lagtime == 0) {
            // Ready to propagate CTMSM, update state and diffusion coefficients
            part.propagateTMSM = true;
            part.setState(part.nextState);
            part.setD(tmsm.Dlist[part.state]);
            part.setDrot(tmsm.Drotlist[part.state]);
        } else {
            part.propagateTMSM = false;
        };
    };
};


/* Integrates rotation/translation and Markovian switch of one particle, with pair interactions
 * (visible only inside the class) */
template<>
void odLangevinMarkovSwitch<ctmsm>::integrateOneMS(int partIndex, std::vector<particleMS> &parts, double timestep) {
    auto &part = parts[partIndex];
    // Do diffusion/rotation propagation taking MSM/CTMSM into account
    double resdt;
    // propagate CTMSM when synchronized and update diffusion coefficients
    part.tcount = 0;
    // Runs one timestep dt, if lagtime < dt, propagates CTMSM as many times as needed
    if (part.lagtime <= timestep) {
        // Run loop until integration for one timestep dt is achieved
        while (part.tcount < timestep) {
            // Propagates MSM only when diffusion and rotation are in sync
            if (part.propagateTMSM) {
                tmsm.propagateNoUpdate(part, 1); // Diffusion coefficients don't need to be updated until dt reaches lagtime.
            }
            // Integrates for one lagtime as long as integration is still under dt
            if ( part.tcount + part.lagtime < timestep ) {
                integrateOne(partIndex, parts, part.lagtime);
                part.tcount += part.lagtime;
                part.setLagtime(0);
                // Ready to propagate MSM, update state and diffusion coefficients
                part.propagateTMSM = true;
                part.setState(part.nextState);
                part.setD(tmsm.Dlist[part.state]);
                part.setDrot(tmsm.Drotlist[part.state]);
                // If current lagtime overtakes dt, integrate up to dt (by resdt) and reset lagtime to remaining portion
            } else {
                resdt = timestep - part.tcount;
                integrateOne(partIndex, parts, resdt);
                part.setLagtime(part.lagtime + part.tcount - timestep);
                part.tcount += resdt; // this means part.tcount = dt, so will exit while loop.
                // If lag time = 0 CTMSM must propagate in next step, otherwise it needs to integrate remaining lagtime.
                if (part.lagtime == 0) {
                    // Ready to propagate CTMSM, update state and diffusion coefficients
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
        integrateOne(partIndex, parts, timestep);
        part.setLagtime(part.lagtime - timestep);
        // If lag time = 0 MSM must propagate in next step, otherwise it needs to integrate remaining lagtime.
        if (part.lagtime == 0) {
            // Ready to propagate CTMSM, update state and diffusion coefficients
            part.propagateTMSM = true;
            part.setState(part.nextState);
            part.setD(tmsm.Dlist[part.state]);
            part.setDrot(tmsm.Drotlist[part.state]);
        } else {
            part.propagateTMSM = false;
        };
    };
};


/*
 * Next two functions should remain at end of file to avoid instantiation before specialization
 */

// Same as integrateOne, but updates time and publicly visible
template<>
void odLangevinMarkovSwitch<ctmsm>::integrate(particleMS &part) {
    integrateOneMS(part, dt);
    clock += dt;
};

// Integrates list of particlesMS (needs to override parent function because it is template based and uses particleMS)
template<>
void odLangevinMarkovSwitch<ctmsm>::integrateList(std::vector<particleMS> &parts, bool pairInteractions) {
    if (pairInteractions) {
        for (int i = 0; i < parts.size(); i++) {
            // Calls version of integrateOne for pair interactions (needs whole particle list)
            integrateOne(i, parts, dt);
        }
    }
    else {
        for (int i = 0; i < parts.size(); i++) {
            // Calls version of integrateOne without pair interactions
            integrateOne(parts[i], dt);
        }
    }
    clock += dt;
};


