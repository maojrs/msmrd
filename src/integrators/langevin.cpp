//
// Created by maojrs on 4/19/21.
//

#include "integrators/langevin.hpp"

namespace msmrd {
    /**
     * Implementation of Lanegvin dynamics integrator class
     * @param dt time step
     * @param seed random generator seed (Note seed = -1 corresponds to random device)
     * @param rotation boolean to indicate if rotational degrees of freedom should be integrated
     */
    langevin::langevin(double dt, long seed, std::string particlesbodytype, std::string integratorScheme)
            : integrator(dt, seed, particlesbodytype), integratorScheme(integratorScheme) {
        rotation = false;
        velocityIntegration = true;
        if (particlesbodytype != "point") {
            throw std::invalid_argument("Langevin integrator only implemented for point particles "
                                        "(no rotation allowed)");
        }
    }

    langevin::langevin(double dt, long seed, std::string particlesbodytype)
            : integrator(dt, seed, particlesbodytype) {
        rotation = false;
        velocityIntegration = true;
        integratorScheme = "BAOAB";
        if (particlesbodytype != "point") {
            throw std::invalid_argument("Langevin integrator only implemented for point particles "
                                        "(no rotation allowed)");
        }
    }

    // Integrate list of particles (need to override in case of MS particles)
    void langevin::integrate(std::vector<particle> &parts) {
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
        if (integratorScheme == "ABOBA") {
            integrateABOBA(parts, dt);
        } else {
            integrateBAOAB(parts, dt);
        }

        // Updates time
        clock += dt;
    }

    // Integrates velocity half a time step given potential or force term
    void langevin::integrateB(std::vector<particle> &parts, double deltat) {
        vec3<double> force;
        for (int i = 0; i < parts.size(); i++) {
            force = 1.0 * forceField[i];
            auto mass = parts[i].mass;
            auto newVel = parts[i].nextVelocity + deltat * force / mass;
            parts[i].setNextVelocity(newVel);
        }
    }

    // Integrates position half a time step given velocity term
    void langevin::integrateA(std::vector<particle> &parts, double deltat) {
        for (int i = 0; i < parts.size(); i++) {
            auto newPos = parts[i].nextPosition + deltat * parts[i].nextVelocity;
            parts[i].setNextPosition(newPos);
        }
    }

    // Integrates velocity full time step given friction and noise term
    void langevin::integrateO(std::vector<particle> &parts, double deltat) {
        for (int i = 0; i < parts.size(); i++) {
            auto eta = KbTemp / parts[i].D; // friction coefficient
            auto mass = parts[i].mass;
            auto xi = std::sqrt(KbTemp * mass * (1 - std::exp(-2 * eta * deltat / mass))) / mass;
            auto newVel = std::exp(-deltat * eta / mass) * parts[i].nextVelocity + xi * randg.normal3D(0, 1);
            parts[i].setNextVelocity(newVel);
        }
    }

    void langevin::integrateBAOAB(std::vector<particle> &parts, double timestep) {
        /* Integrate BAOA part and save next positions/orientations in parts[i].next. Note final
         * position already obtained at the end of these integrations. */
        integrateB(parts, dt/2.0);
        integrateA(parts, dt/2.0);
        integrateO(parts, dt);
        integrateA(parts, dt/2.0);

        // Enforce boundary; sets new positions into parts[i].nextPosition (only if particle is active)
        enforceBoundary(parts);

        // Update particlesposition to recalculate force for last scheme step "B"
        updatePositionOrientation(parts);

        // Recalculate force field (to be used also in next time step)
        calculateForceTorqueFields(parts);

        // Integrate final B part with newly calculated force to obtain final velocities
        integrateB(parts, dt/2.0);

        // Update velocity based on parts[i].nextVelocity
        updateVelocities(parts);
    };

    void langevin::integrateABOBA(std::vector<particle> &parts, double timestep) {

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