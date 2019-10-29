//
// Created by maojrs on 8/16/18.
//

#include "integrators/overdampedLangevinMarkovSwitch.hpp"

namespace msmrd {
   /**
    * Implementation of over-damped Langevin dynamics with Markovian switch integrator class
    * (constructors defined in headers since templates were used).
    */

   // Aliases for classes with long names
    using msm = msmrd::discreteTimeMarkovStateModel;
    using ctmsm = msmrd::continuousTimeMarkovStateModel;


    /* Integrates diffusion and rotation of one particle, called by the integrateOneMS (visible only
     * inside the class). Note it cannot be inherited from overdampedLangevin::integrateOne since it
     * uses a particle list instead of a particle list. Using list of pointers possible but not ideal
     * due to dependencies. */
    template<>
    void overdampedLangevinMarkovSwitch<ctmsm>::integrateOne(int partIndex, std::vector<particle> &parts, double timestep) {
        vec3<double> force;
        vec3<double> torque;
        force = forceField[partIndex];
        torque = torqueField[partIndex];
        translate(parts[partIndex], force, timestep);
        if (rotation) {
            rotate(parts[partIndex], torque, timestep);
        }
    };


    /* Integrates rotation/translation and Markovian switch of one particle, with pair interactions
     * (visible only inside the class) */
    template<>
    void overdampedLangevinMarkovSwitch<ctmsm>::integrateOneMS(int partIndex, std::vector<particle> &parts, double timestep) {
        auto &part = parts[partIndex];
        auto tmsm = MSMlist[part.type];
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
                if (part.tcount + part.lagtime < timestep) {
                    integrateOne(partIndex, parts, part.lagtime);
                    part.tcount += part.lagtime;
                    part.setLagtime(0);
                    // Ready to propagate MSM, update state and diffusion coefficients
                    part.propagateTMSM = true;
                    part.updateState(); // Sets calculated next state as current state
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
                        part.updateState(); // Sets calculated next state as current state
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
                part.updateState(); // Sets calculated next state as current state
                part.setD(tmsm.Dlist[part.state]);
                part.setDrot(tmsm.Drotlist[part.state]);
            } else {
                part.propagateTMSM = false;
            };
        };
    };


    /*
     * Next function should remain at end of file to avoid instantiation before specialization
     */

    /* Integrates list of particle particles (needs to override parent function because it is
     * template based and uses particle) */
    template<>
    void overdampedLangevinMarkovSwitch<ctmsm>::integrate(std::vector<particle> &parts) {
        // Calculate forces and torques and save them into forceField and torqueField
        calculateForceTorqueFields<particle>(parts);

        // Integrate and save next positions/orientations in parts[i].next***
        for (int i = 0; i < parts.size(); i++) {
            if (parts[i].activeMSM) {
                integrateOneMS(i, parts, dt);
            } else {
                integrateOne(i, parts, dt);
            }
        }
        // Enforce boundary and set new positions into parts[i].nextPosition
        enforceBoundary(parts);

        /* Update positions and orientations (sets calculated next position/orientation
         * calculated by integrator and boundary as current position/orientation). */
        updatePositionOrientation(parts);

        clock += dt;
    }

}


