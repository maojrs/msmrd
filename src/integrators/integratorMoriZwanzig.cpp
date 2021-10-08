//
// Created by maojrs on 6/14/21.
//

#include "integrators/integratorMoriZwanzig.hpp"

namespace msmrd {

    /* Loads auxiliary variable into raux. In this case this correspond to the potential and noise term,
     * which can be calculated as M(dV) + gamma V_i dt .*/
    void integratorMoriZwanzig::loadAuxiliaryValues(std::vector<particle> &parts) {
        for (int i = 0; i < parts.size(); i++) {
            // If particle type corresponds to one of the distinguished particles, sample its value
            if (std::find(distinguishedTypes.begin(), distinguishedTypes.end(),
                          parts[i].type) != distinguishedTypes.end()) {
                double eta = KbTemp / parts[i].D;
                auto momentumDiff = parts[i].mass * (parts[i].nextVelocity - parts[i].velocity);
                auto raux = momentumDiff + eta * parts[i].velocity * dt;
                parts[i].raux = raux;
            }
        }
    };

    void integratorMoriZwanzig::integrate(std::vector<particle> &parts) {
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

        // Load auxiliary variables into distinguished particle before updating velocities
        loadAuxiliaryValues(parts);

        // Update velocity based on parts[i].nextVelocity
        updateVelocities(parts);

        // Updates time
        clock += dt;
    };

}
