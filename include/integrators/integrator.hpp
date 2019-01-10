//
// Created by dibakma on 27.06.18.
//

#pragma once
#include <array>
#include <utility>
#include <memory>
#include "boundaries/boundary.hpp"
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
        randomgen randg;

        std::vector<vec3<double>> forceField;
        std::vector<vec3<double>> torqueField;

        bool boundaryActive = false;
        bool externalPotentialActive = false;
        bool pairPotentialActive = false;

        // Boundary pointer
        boundary *domainBoundary;

        // External potentials pointers (externalPot)
        externalPotential *externalPot;

        // Pair potentials pointers (pairPot)
        pairPotential *pairPot;


        /**
        * @param KbTemp = Boltzman constant times temperature
        * @param dt time step
        * @param seed variable for random number generation (Note seed = -1 corresponds to random device)
        * @param particlesbodytype body type of particles to integrate. It determines rotation integrator behavior, can
        * be either point, rod or rigidbody, and it is determined by orientational degrees of freedom, points
        * have no orientation, rods need only one vector and rigidsolid requires a complete quaternion). The particles
        * body type can also be pointMix, rodMix or rigidbodyMix, where there are different types of
        * particles (particleMS.type) with the same body type. This allows the integrator to pass the particle type to
        * the potential function. Note the particlesbodytype has to match the potential used.
        * @param rotation boolean to indicate if rotation should be integrated
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
        */

        // Protected abstract functions required by main integrate public routine
        virtual void integrateOne(int index, std::vector<particle> &parts, double timestep) = 0;

        virtual void translate(particle &part, vec3<double> force, double dt) = 0;

        virtual void rotate(particle &part, vec3<double> torque, double dt) = 0;


        /* Protected functions to get forces and torques due to external or pair potentials for integrator.
         * The template PARTICLE can take values of particle or particleMS or other custom defined particles.
         * This is required because vectors of child classes are not recognized as childs of vector of parents class.*/
        template< typename PARTICLE >
        void calculateForceTorqueFields(const std::vector<PARTICLE> &parts); // external + pairs

        template< typename PARTICLE >
        void calculateExternalForceTorques(const std::vector<PARTICLE> &parts, int numParticles);

        template< typename PARTICLE >
        void calculatePairsForceTorques(const std::vector<PARTICLE> &parts, int numParticles);


    public:
        double clock;
        /**
         * @param clock keeps track of global time
         */

        integrator(double dt, long seed, std::string particlesbodytype);

        // Main public functions definitions
        void integrate(std::vector<particle> &parts);

        double getClock() const { return clock; }

        void setKbT(double kbt) { KbTemp = kbt; }


        // Potential and boundary related functions
        void setBoundary(boundary *bndry);

//        void setExternalPotential(externalPotential<particle> *pot);
//        void setExternalPotentialMS(externalPotential<particleMS> *pot);
//        void setPairPotential(externalPotential<particle> *pot);
//        void setPairPotentialMS(externalPotential<particleMS> *pot);


        void setExternalPotential(externalPotential *pot);

        void setPairPotential(pairPotential *pot);

    };


    /**
     * Implements template functions to calculate force and torque fields
     * Note particle can be the particle class or any of its child classes, like particleMS
     * @param parts is a vector of PARTICLE objects
     */

    /* Calculates the total force and torque fields experienced by all particles, including external and
     * interaction pairs potentials. Relies on two other auxiliary functions to divide the work. */
    template <typename PARTICLE>
    void integrator::calculateForceTorqueFields(const std::vector<PARTICLE> &parts) {
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


    /* Calculate external forces and torques due to interactiion with an external potential and save it into
     * forceField and torqueField. Note different bodytypes require different fucntion calls */
    template <typename PARTICLE>
    void integrator::calculateExternalForceTorques(const std::vector<PARTICLE> &parts, int numParticles) {
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
    void integrator::calculatePairsForceTorques(const std::vector<PARTICLE> &parts, int numParticles) {
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

}




