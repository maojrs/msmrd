//
// Created by maojrs on 4/19/21.
//

#include "integrators/integratorLangevin.hpp"

namespace msmrd {
    /**
     * Implementation of Lanegvin dynamics integrator virtual class
     * @param dt time step
     * @param seed random generator seed (Note seed = -1 corresponds to random device)
     * @param rotation boolean to indicate if rotational degrees of freedom should be integrated
     * @param frictionCoefficient friction coefficient in Langevin equation in units of mass/time.
     */

    integratorLangevin::integratorLangevin(double dt, long seed, std::string particlesbodytype, double Gamma)
            : integrator(dt, seed, particlesbodytype), Gamma(Gamma) {
        rotation = false;
        if (particlesbodytype != "point") {
            throw std::invalid_argument("Langevin integrator only implemented for point particles "
                                        "(no rotation allowed)");
        }
    }

    // Integrate list of particles (need to override in case of MS particles)
    void integratorLangevin::integrate(std::vector<particle> &parts) {
        vec3<double> force;
        vec3<double> torque;

        // Calculate forces and torques and save them into forceField and torqueField for the first run
        if (firstRun) {
            calculateForceTorqueFields(parts);
            firstRun = false;
        }

        /* As integrateB,A and O functions are used recursively, they act on the nextPosition and nextVelocity
        * variable, so we need to set nextVariables equal to current ones. */
        for (int i = 0; i < parts.size(); i++) {
            parts[i].setNextPosition(parts[i].position);
            parts[i].setNextVelocity(parts[i].velocity);
        }

        // Integrates one time step of the chosen integrator, default BAOAB
        integrateOneTimestep(parts, dt);

        // Updates time
        clock += dt;
    }

    // Integrates velocityfor deltat given potential or force term
    void integratorLangevin::integrateB(std::vector<particle> &parts, double deltat) {
        vec3<double> force;
        for (int i = 0; i < parts.size(); i++) {
            force = 1.0 * forceField[i];
            auto mass = parts[i].mass;
            auto newVel = parts[i].nextVelocity + deltat * force / mass;
            parts[i].setNextVelocity(newVel);
        }
    }

    // Integrates position for deltat time step given velocity term
    void integratorLangevin::integrateA(std::vector<particle> &parts, double deltat) {
        for (int i = 0; i < parts.size(); i++) {
            auto newPos = parts[i].nextPosition + deltat * parts[i].nextVelocity;
            parts[i].setNextPosition(newPos);
        }
    }

    // Integrates velocity for timstep delta t given friction and noise term
    void integratorLangevin::integrateO(std::vector<particle> &parts, double deltat) {
        for (int i = 0; i < parts.size(); i++) {
            auto mass = parts[i].mass;
            auto xi = std::sqrt(KbTemp * mass * (1 - std::exp(-2 * Gamma * deltat / mass))) / mass;
            auto newVel = std::exp(-deltat * Gamma / mass) * parts[i].nextVelocity + xi * randg.normal3D(0, 1);
            parts[i].setNextVelocity(newVel);
        }
    }


    /*
     * Implementation of different Langevin integrators using different schemes, eahc in a different class
     */

    void langevinSemiImplicitEuler::integrateOneTimestep(std::vector<particle> &parts, double deltat) {
        // Calculate force field
        calculateForceTorqueFields(parts);

        /* Integrate OBA part and save next positions/orientations in parts[i].next. Corresponds to semi-implicit
         * Euler method. */
        integrateO(parts, deltat);
        integrateB(parts, deltat);
        integrateA(parts, deltat);

        // Enforce boundary; sets new positions into parts[i].nextPosition (only if particle is active)
        enforceBoundary(parts);

        // Update particlesposition to recalculate force for last scheme step "B"
        updatePositionOrientation(parts);

        // Update velocity based on parts[i].nextVelocity
        updateVelocities(parts);
    };

    void langevinBAOAB::integrateOneTimestep(std::vector<particle> &parts, double deltat) {
        /* Integrate BAOA part and save next positions/orientations in parts[i].next. Note final
         * position already obtained at the end of these integrations. */
        integrateB(parts, deltat/2.0);
        integrateA(parts, deltat/2.0);
        integrateO(parts, deltat);
        integrateA(parts, deltat/2.0);

        // Enforce boundary; sets new positions into parts[i].nextPosition (only if particle is active)
        enforceBoundary(parts);

        // Update particlesposition to recalculate force for last scheme step "B"
        updatePositionOrientation(parts);

        // Recalculate force field (to be used also in next time step)
        calculateForceTorqueFields(parts);

        // Integrate final B part with newly calculated force to obtain final velocities
        integrateB(parts, deltat/2.0);

        // Update velocity based on parts[i].nextVelocity
        updateVelocities(parts);
    };

    void langevinABOBA::integrateOneTimestep(std::vector<particle> &parts, double timestep) {

        integrateA(parts, dt/2.0);

        // Update particlesposition to recalculate force for next step "B"
        updatePositionOrientation(parts);

        // Calculate force field
        calculateForceTorqueFields(parts);

        integrateB(parts, dt/2.0);
        integrateO(parts, dt);
        integrateB(parts, dt/2.0);
        integrateA(parts, dt/2.0);

        // Enforce boundary; sets new positions into parts[i].nextPosition (only if particle is active)
        enforceBoundary(parts);

        // Update particlesposition to recalculate force for last scheme step "B"
        updatePositionOrientation(parts);

        // Update velocity based on parts[i].nextVelocity
        updateVelocities(parts);
    };





}