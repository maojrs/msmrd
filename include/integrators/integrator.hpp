//
// Created by dibakma on 27.06.18.
//

#pragma once
#include <array>
#include <utility>
#include <memory>
#include "boundaries/boundary.hpp"
#include "boundaries/noBoundary.hpp"
#include "particle.hpp"
#include "randomgen.hpp"
#include "potentials/potentials.hpp"


namespace msmrd {
    /**
     * Abstract integrators base class declaration
     */
    class integrator {
    protected:
        double KbTemp = 1.0;
        double dt;
        long seed;
        std::string particlesbodytype;
        bool rotation = false;
        bool velocityIntegration = false;
        randomgen randg = randomgen();


        std::vector<vec3<double>> forceField;
        std::vector<vec3<double>> torqueField;

        bool boundaryActive = false;
        bool externalPotentialActive = false;
        bool pairPotentialActive = false;

        // Boundary pointer (set default boundary (inactive), i.e. no boundary)
        boundary *domainBoundary;

        // External potentials pointers (externalPot)
        externalPotential *externalPot;

        // Pair potentials pointers (pairPot)
        pairPotential *pairPot;


        /**
        * @param KbTemp = Boltzman constant times temperature (default value assumed to be 1 to use
        * reduced potential and force. Can be adjusted if neccesary.
        * @param dt time step
        * @param seed variable for random number generation (Note seed = -1 corresponds to random device)
        * @param particlesbodytype body type of particles to integrate. It determines rotation integrator behavior, can
        * be either point, rod or rigidbody, and it is determined by orientational degrees of freedom, points
        * have no orientation, rods need only one vector and rigidsolid requires a complete quaternion).
        * @param rotation boolean to indicate if rotation should be integrated
        * @param velocityIntegration boolean to indicate if velocity should be integrated. Can only be set to
        * true for Langevin type integrators.
        * @param randg random number generator based in mt19937
        *
        * @param forceField stores force experienced by each particle at a given time
        * @param torqueField stores torque experienced by each particle at a given time
        *
        * @param boundaryActive indicates if a boundary conditions is active
        * @param externalPotActive indicates if external potential has been set
        * @param pairPotActive indicates if potential potential has been set

        * @param *domainBoundary pointer to the boundary object to be used.
        * Note all following potentials default to zero and not every integrator will make use of all this potentials
        * @param *externalPot pointer to external potential
        * @param *pairPot pointer to pair potential between two particles
        * @param clock keeps track of global time
        */

        // Protected abstract functions required by main integrate public routine
        virtual void integrateOne(int index, std::vector<particle> &parts, double timestep) = 0;

        virtual void translate(particle &part, vec3<double> force, double dt) = 0;

        virtual void rotate(particle &part, vec3<double> torque, double dt) = 0;


        /* Protected functions to get forces and torques due to external or pair potentials for integrator.
         * The template PARTICLE can take values of particle or particleMS or other custom defined particles.
         * This is required because vectors of child classes are not recognized as childs of vector of parents class.*/
        template< typename PARTICLE >
        void calculateForceTorqueFields(std::vector<PARTICLE> &parts); // external + pairs

        template< typename PARTICLE >
        void calculateExternalForceTorques(std::vector<PARTICLE> &parts, int numParticles);

        template< typename PARTICLE >
        void calculatePairsForceTorques(std::vector<PARTICLE> &parts, int numParticles);

        // Other functions used by most integrators, so defined here as template functions

        template< typename PARTICLE >
        void updatePositionOrientation(std::vector<PARTICLE> &parts);

        template< typename PARTICLE >
        void updateVelocities(std::vector<PARTICLE> &parts);

        template< typename PARTICLE >
        void enforceBoundary(std::vector<PARTICLE> &parts);


    public:
        double clock;

        integrator(double dt, long seed, std::string particlesbodytype);

        // Main public functions definitions
        virtual void integrate(std::vector<particle> &parts);

        vec3<double> calculateRelativePosition(vec3<double> p1, vec3<double> p2);


        // Getters and setters
        void setBoundary(boundary *bndry);

        void setExternalPotential(externalPotential *pot);

        void setPairPotential(pairPotential *pot);

        void setKbT(double kbt) { KbTemp = kbt; }

        double getClock() const { return clock; }

        void setClock(double newTime) { clock = newTime; }

        void resetClock() {clock = 0.0;}

        std::string getParticlesBodyType() const { return particlesbodytype; }

        bool isBoundaryActive() { return boundaryActive; }

        boundary* getBoundary() { return domainBoundary; }

    };


    /**
     * Implements template functions to calculate force and torque fields
     * Note particle can be the particle class or any of its child classes, like particleMS
     * @param parts is a vector of PARTICLE objects
     */

    /* Calculates the total force and torque fields experienced by all particles, including external and
     * interaction pairs potentials. Relies on two other auxiliary functions to divide the work. */
    template <typename PARTICLE>
    void integrator::calculateForceTorqueFields(std::vector<PARTICLE> &parts) {
        unsigned int N = static_cast<int>(parts.size());

        // Resize force/torque fields array if the number of particles have changed
        if (static_cast<int>(forceField.size()) != N) {
            forceField.resize(N);
            torqueField.resize(N);
        }

        /* Add forces and torques due to external potential from all particles to fields,
         * initialize to zero if there is no external potential. */
        if (externalPotentialActive) {
            calculateExternalForceTorques(parts, N);
        } else {
            for (int i = 0; i < N; i++) {
                forceField[i] = vec3<double> (0,0,0);
                torqueField[i] = vec3<double> (0,0,0);
            }
        }

        // Add forces and torques coming from pair interactions.
        if (pairPotentialActive) {
            calculatePairsForceTorques(parts, N);
        }
    }


    /* Calculate external forces and torques due to interaction with an external potential and save it into
     * forceField and torqueField. Note different bodytypes require different fucntion calls */
    template <typename PARTICLE>
    void integrator::calculateExternalForceTorques(std::vector<PARTICLE> &parts, int numParticles) {
        std::array<vec3<double>, 2> forctorq;
        for (int i = 0; i < numParticles; i++) {
            forctorq = externalPot->forceTorque(parts[i]);
            forceField[i] = 1.0 * forctorq[0];
            torqueField[i] = 1.0 * forctorq[1];
        }
    };

    /* Calculate pairs forces and torques due to interaction with other particles and sums it into
     * forceField and torqueField. Note different bodytypes require different function calls. Also note
     * to avoid duplicate calls to pairPotential.forceTorque, the two loops cover all the pair interactions
     * only once. The function pairPotential.forceTorque returns the force and torque exerted on particle1
     * and the force and torque exerted on particle2 (in that order), from their mutual interaction. */
    template <typename PARTICLE>
    void integrator::calculatePairsForceTorques(std::vector<PARTICLE> &parts, int numParticles) {
        unsigned int N = static_cast<int>(parts.size());
        std::array<vec3<double>, 4> forctorq;
        // Calculate the forces and torque for each possible interaction
        for (int i = 0; i < numParticles; i++) {
            for (int j = i + 1; j < numParticles; j++) {
                forctorq = pairPot->forceTorque(parts[i], parts[j]);
                forceField[i] += 1.0*forctorq[0];
                torqueField[i] += 1.0*forctorq[1];
                forceField[j] += 1.0*forctorq[2];
                torqueField[j] += 1.0*forctorq[3];
            }
        }
    }


    /* Update positions and orientations (sets calculated next position/orientation
     * calculated by integrator and boundary as current position/orientation). Only
     * updated if particle is active. Orientation only updated if rotation is active */
    template <typename PARTICLE>
    void integrator::updatePositionOrientation(std::vector<PARTICLE> &parts){
        for (int i = 0; i < parts.size(); i++) {
            if (parts[i].isActive()) {
                parts[i].updatePosition();
                if (rotation) {
                    parts[i].updateOrientation();
                }
            }
        }
    }

    /* Update velocities (sets calculated next velocity
     * calculated by integrator and boundary as current velocity). Only
     * updated if particle is active and if velocityIntegration is active. */
    template <typename PARTICLE>
    void integrator::updateVelocities(std::vector<PARTICLE> &parts) {
        if (velocityIntegration) {
            for (int i = 0; i < parts.size(); i++) {
                if (parts[i].isActive()) {
                    parts[i].updateVelocity();
                }
            }
        }
    }

    /* Enforces boundary; sets new positions given by boundary conditions (e.g. periodic bounary)
     * into parts[i].nextPosition. This is only done if particle is active (which is the default
     * value for a newly created particle) */
    template< typename PARTICLE >
    void integrator::enforceBoundary(std::vector<PARTICLE> &parts) {
        for (auto &part : parts) {
            if (part.isActive() and boundaryActive) {
                domainBoundary->enforceBoundary(part);
            }
        }
    };

}




