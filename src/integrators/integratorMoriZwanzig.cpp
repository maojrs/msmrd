//
// Created by maojrs on 6/14/21.
//

#include "integrators/integratorMoriZwanzig.hpp"

namespace msmrd {

    integratorMoriZwanzig::integratorMoriZwanzig(double dt, long seed, std::string particlesbodytype, double Gamma) :
    langevinABOBA(dt, seed, particlesbodytype, Gamma) {};

    /* Loads auxiliary variable into raux. In this case this correspond to the potential and noise term,
     * which can be calculated as M(dV) + Gamma V_i dt .*/
    void integratorMoriZwanzig::loadAuxiliaryValues(std::vector<particle> &parts,
                                                    std::vector<vec3<double>> auxForces) {
        for (int i = 0; i < parts.size(); i++) {
            // If particle type corresponds to one of the distinguished particles, sample its value
            if (std::find(distinguishedTypes.begin(), distinguishedTypes.end(),
                          parts[i].type) != distinguishedTypes.end()) {
                auto mass = parts[i].mass;
                auto interactionTerm = auxForces[i]*dt/(2*mass) * (1 + std::exp(-dt * Gamma / mass));
                auto noiseTerm = parts[i].raux2;
                parts[i].raux = interactionTerm + noiseTerm;
            }
        }
    };

    void integratorMoriZwanzig::integrateOneTimestep(std::vector<particle> &parts, double timestep) {

        integrateA(parts, dt/2.0);

        // Update particlesposition to recalculate force for next step "B"
        updatePositionOrientation(parts);

        // Calculate force field
        calculateForceTorqueFields(parts);

        integrateB(parts, dt/2.0);
        integrateO(parts, dt);
        integrateB(parts, dt/2.0);
        integrateA(parts, dt/2.0);

        // Load auxiliary variables into distinguished particle before updating velocities
        loadAuxiliaryValues(parts, auxForceField);

        // Enforce boundary; sets new positions into parts[i].nextPosition (only if particle is active)
        enforceBoundary(parts);

        // Update particlesposition to recalculate force for last scheme step "B"
        updatePositionOrientation(parts);

        // Update velocity based on parts[i].nextVelocity
        updateVelocities(parts);
    };

    /* Integrates velocity full time step given friction and noise term, saves noise term in raux2 variable.
     * Specialized version of the one implemented in the integratorLangevin integrator. */
    void integratorMoriZwanzig::integrateO(std::vector<particle> &parts, double deltat) {
        for (int i = 0; i < parts.size(); i++) {
            auto mass = parts[i].mass;
            auto xi = std::sqrt((KbTemp/mass) * (1 - std::exp(-2 * Gamma * deltat / mass)));
            auto noiseTerm = xi * randg.normal3D(0, 1);
            auto newVel = std::exp(-deltat * Gamma / mass) * parts[i].nextVelocity + noiseTerm;
            parts[i].setNextVelocity(newVel);
            parts[i].raux2 = noiseTerm;
        }
    }


    /*
     * Class implementations for integratorMoriZwanzig2, alternative version.
     */

    /* Loads auxiliary variable into raux. In this case this correspond to the potential and noise term,
     * which can be calculated as M(dV) + Gamma V_i dt .*/
    void integratorMoriZwanzig2::loadAuxiliaryValues(std::vector<particle> &parts,
                                                    std::vector<vec3<double>> auxForces) {
        for (int i = 0; i < parts.size(); i++) {
            // If particle type corresponds to one of the distinguished particles, sample its value
            if (std::find(distinguishedTypes.begin(), distinguishedTypes.end(),
                          parts[i].type) != distinguishedTypes.end()) {
                auto mass = parts[i].mass;
                auto interactionTerm = auxForces[i]*dt/(2*mass) * (1 + std::exp(-dt * Gamma / mass));
                parts[i].raux = interactionTerm;
                // Noise term (parts[i].raux2) loaded on integrateO.
            }
        }
    };


    /*
     * Class implementations for integratorMoriZwanzig2, alternative version.
     */

    /* Loads auxiliary variable into raux as in original function, but only saves the index-component. */
    void integratorMoriZwanzigConstrained1D::loadAuxiliaryValues(std::vector<particle> &parts,
                                                    std::vector<vec3<double>> auxForces) {
        for (int i = 0; i < parts.size(); i++) {
            // If particle type corresponds to one of the distinguished particles, sample its value
            if (std::find(distinguishedTypes.begin(), distinguishedTypes.end(),
                          parts[i].type) != distinguishedTypes.end()) {
                auto mass = parts[i].mass;
                auto interactionTerm = auxForces[i]*dt/(2*mass) * (1 + std::exp(-dt * Gamma / mass));
                auto noiseTerm = parts[i].raux2;
                parts[i].raux[0] = interactionTerm[0] + noiseTerm[0];
                parts[i].raux[1] = 0.0;
                parts[i].raux[2] = 0.0;
            }
        }
    };

    /* Integrates Mori Zwanzig constraining the dimer in the x-axis
     */
    void integratorMoriZwanzigConstrained1D::integrateOneTimestep(std::vector<particle> &parts, double timestep) {

        integrateA(parts, dt/2.0);

        // Update particlesposition to recalculate force for next step "B"
        updatePositionOrientation(parts);

        // Calculate force field
        calculateForceTorqueFields(parts);

        integrateB(parts, dt/2.0);
        integrateO(parts, dt);
        integrateB(parts, dt/2.0);
        integrateA(parts, dt/2.0);

        // Load auxiliary variables into distinguished particle before updating velocities
        loadAuxiliaryValues(parts, auxForceField);

        // Enforce boundary; sets new positions into parts[i].nextPosition (only if particle is active)
        enforceBoundary(parts);

        // Update particles position to recalculate force for last scheme step "B", only in one dimension for the dimer
        updatePositionOrientation1D(parts);

        // Update velocity based on parts[i].nextVelocity, only in one dimension for the dimer
        updateVelocities1D(parts);
    };

    void integratorMoriZwanzigConstrained1D::updatePositionOrientation1D(std::vector<particle> &parts){
        for (int i = 0; i < parts.size(); i++) {
            if (parts[i].isActive()) {
                if (std::find(distinguishedTypes.begin(), distinguishedTypes.end(),
                              parts[i].type) != distinguishedTypes.end()) {
                    // Update only te component index of position and orientation
                    parts[i].updatePositionIndex(index1D);
                    if (rotation) {
                        parts[i].updateOrientationIndex(index1D);
                    }
                } else {
                    parts[i].updatePosition();
                    if (rotation) {
                        parts[i].updateOrientation();
                    }
                }
            }
        }
    }

    void integratorMoriZwanzigConstrained1D::updateVelocities1D(std::vector<particle> &parts) {
        for (int i = 0; i < parts.size(); i++) {
            if (parts[i].isActive()) {
                if (std::find(distinguishedTypes.begin(), distinguishedTypes.end(),
                              parts[i].type) != distinguishedTypes.end()) {
                    parts[i].updateVelocityIndex(index1D);
                } else {
                    parts[i].updateVelocity();
                }
            }
        }
    }


}
