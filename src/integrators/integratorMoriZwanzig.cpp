//
// Created by maojrs on 6/14/21.
//

#include "integrators/integratorMoriZwanzig.hpp"

namespace msmrd {

    integratorMoriZwanzig::integratorMoriZwanzig(double dt, long seed, std::string particlesbodytype, double frictionCoefficient) :
    langevin(dt, seed, particlesbodytype, frictionCoefficient, "modifiedABOBA") {};

    /* Loads auxiliary variable into raux. In this case this correspond to the potential and noise term,
     * which can be calculated as M(dV) + gamma V_i dt .*/
    void integratorMoriZwanzig::loadAuxiliaryValues(std::vector<particle> &parts,
                                                    std::vector<vec3<double>> pairsForces) {
        for (int i = 0; i < parts.size(); i++) {
            // If particle type corresponds to one of the distinguished particles, sample its value
            if (std::find(distinguishedTypes.begin(), distinguishedTypes.end(),
                          parts[i].type) != distinguishedTypes.end()) {
                auto eta = KbTemp / parts[i].D;
                auto mass = parts[i].mass;
                auto interactionTerm = pairsForces[i]*dt/(2*mass) * (1 + std::exp(-dt * eta / mass));
                auto noiseTerm = parts[i].raux2;
                parts[i].raux = interactionTerm + noiseTerm;
            }
        }
    };

    void integratorMoriZwanzig::integrate(std::vector<particle> &parts) {
        vec3<double> force;
        vec3<double> torque;
        unsigned int N = static_cast<int>(parts.size());

        // Calculate forces and torques and save them into forceField and torqueField for the first run
        if (firstRun) {
            calculateForceTorqueFields(parts);
            firstRun = false;
        }

        /* As integrate B,A and O functions are used recursively they act on the nextPosition and nextVelocity
         * variable, so we need to set nextVariables equal to current ones. */
        for (int i = 0; i < parts.size(); i++) {
            parts[i].setNextPosition(parts[i].position);
            parts[i].setNextVelocity(parts[i].velocity);
        }

        integrateA(parts, dt/2.0);

        // Enforce boundary; sets new positions into parts[i].nextPosition (only if particle is active)
        enforceBoundary(parts);

        // Update particlesposition to recalculate force for next step "B"
        updatePositionOrientation(parts);

        /* Explicit calculation of force torque fields. Required to extract value of external potential
         * for raux variables. */
        if (externalPotentialActive) {
            calculateExternalForceTorques(parts, N);
        }

        // Save external potential in temporary variable
        auto externalForce = forceField;

        if (pairPotentialActive) {
            calculatePairsForceTorques(parts, N);
        }

        // Save external potential in temporary variable
        std::vector<vec3<double>>  pairsForces;
        pairsForces.resize(N);
        for (int i = 0; i < pairsForces.size(); i++) {
            pairsForces[i] = forceField[i] - externalForce[i];
        }

        integrateB(parts, dt/2.0);
        integrateO(parts, dt);
        integrateB(parts, dt/2.0);
        integrateA(parts, dt/2.0);

        // Load auxiliary variables into distinguished particle before updating velocities
        loadAuxiliaryValues(parts, pairsForces);

        // Enforce boundary; sets new positions into parts[i].nextPosition (only if particle is active)
        enforceBoundary(parts);

        // Update particlesposition to recalculate force for last scheme step "B"
        updatePositionOrientation(parts);

        // Update velocity based on parts[i].nextVelocity
        updateVelocities(parts);

        // Updates time
        clock += dt;
    };


    /* Integrates velocity full time step given friction and noise term, svaes noise term in raux2 variable.
     * Specialized version of the one implemented in the langevin integrator. */
    void integratorMoriZwanzig::integrateO(std::vector<particle> &parts, double deltat) {
        auto eta = frictionCoefficient;
        for (int i = 0; i < parts.size(); i++) {
            auto mass = parts[i].mass;
            auto xi = std::sqrt((KbTemp/mass) * (1 - std::exp(-2 * eta * deltat / mass)));
            auto noiseTerm = xi * randg.normal3D(0, 1);
            auto newVel = std::exp(-deltat * eta / mass) * parts[i].nextVelocity + noiseTerm;
            parts[i].setNextVelocity(newVel);
            parts[i].raux2 = noiseTerm;
        }
    }

}
