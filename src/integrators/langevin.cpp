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
    langevin::langevin(double dt, long seed, std::string particlesbodytype)
            : integrator(dt, seed, particlesbodytype) {
        rotation = false;
        velocityIntegration = true;
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

        /* Integrate BAOA part and save next positions/orientations in parts[i].next. As integrateB,A and O
         * functions are used recursively they act on the nextPosition and nextVelocity variable. Note final
         * position already obtained at the end of this loop. */
        for (int i = 0; i < parts.size(); i++) {
            parts[i].setNextPosition(parts[i].position);
            parts[i].setNextVelocity(parts[i].velocity);
            integrateB(i, parts, dt);
            integrateA(i, parts, dt);
            integrateO(i, parts, dt);
            integrateA(i, parts, dt);
        }

        // Enforce boundary; sets new positions into parts[i].nextPosition (only if particle is active)
        enforceBoundary(parts);

        // Update particlesposition to recalculate force for last scheme step "B"
        updatePositionOrientation(parts);

        // Recalculate force field (to be used also in next time step)
        calculateForceTorqueFields(parts);

        // Integrate final B part with newly calculated force to obtain final velocities
        for (int i = 0; i < parts.size(); i++) {
            integrateB(i, parts, dt);
        }

        // Update velocity based on parts[i].nextVelocity
        updateVelocities(parts);

        // Updates time
        clock += dt;
    }

    // Integrates velocity half a time step given potential or force term
    void langevin::integrateB(int partIndex, std::vector<particle> &parts, double dt) {
        vec3<double> force;
        force = 1.0 * forceField[partIndex];
        auto newVel = parts[partIndex].nextVelocity + dt/2 * force;
        parts[partIndex].setNextVelocity(newVel);
    }

    // Integrates position half a time step given velocity term
    void langevin::integrateA(int partIndex, std::vector<particle> &parts, double dt) {
        auto newPos = parts[partIndex].nextPosition + dt/2 * parts[partIndex].nextVelocity;
        parts[partIndex].setNextPosition(newPos);
    }

    // Integrates velocity full time step given friction and noise term
    void langevin::integrateO(int partIndex, std::vector<particle> &parts, double dt) {
        double eta = KbTemp / parts[partIndex].D; // friction coefficient
        double xi = std::sqrt(KbTemp*(1 - std::exp(-2*eta*dt)));
        auto mass = parts[partIndex].mass;
        auto newVel = std::exp(-dt*eta) * parts[partIndex].nextVelocity + xi/mass * randg.normal3D(0, 1);
        parts[partIndex].setNextVelocity(newVel);
    }



}