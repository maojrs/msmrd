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
        std::string bodytype;
        bool rotation;
        randomgen randg;

        bool boundaryActive = false;
        bool externalPotentialActive = false;
        bool pairPotentialActive = false;

        std::vector<vec3<double>> forceField;
        std::vector<vec3<double>> torqueField;
        std::vector<vec3<double>> forcePairField;
        std::vector<vec3<double>> torquePairField;


        /**
        * @param KbTemp = Boltzman constant times temperature
        * @param dt time step
        * @param seed variable for random number generation (Note seed = -1 corresponds to random device)
        * @param bodytype body type of particles to integrate; determines rotation integrator behavior, can
        * be either point, rod or rigidbody, and it is determined by orientational degrees of freedom, points
        * have no orientation, rods need only one vector and rigidsolid requires a complete quaternion).
        * @param rotation boolean to indicate if rotation should be integrated
        * @param randg random number generator based in mt19937
        * @param boundaryActive indicates if a boundary conditions is active
        * @param externalPotActive indicates if external potential has been set
        * @param pairPotActive indicates if potential potential has been set
        * @param forceField stores forceField for each particle being integrated at the moment
        * @param torqueField stores torqueField for each particle being integrated at the moment
        * @param forcePairField stores all forces between pairs of particles, used to obtain forceField
        * @param torquePairField stores all torques between pairs of particles, used to obtain torqueField
        */

        // Protected abstract functions
        virtual void integrateOne(int index, std::vector<particle> &parts,
                                  double timestep) = 0; // Version for pair interactions (full particle list required)
        virtual void translate(particle &part, vec3<double> force, double dt) = 0;

        virtual void rotate(particle &part, vec3<double> torque, double dt) = 0;

        // Protected functions to get forces and torques due to external or pair potentials for integrator
        template< typename PARTICLE >
        void calculateTotalForceTorqueFields(std::vector<PARTICLE> &parts); // external + pairs

        template< typename PARTICLE >
        void calculatePairsForceTorqueFields(std::vector<PARTICLE> &parts); // only pairs/auxiliary

        // std::array<vec3<double>, 2> getOnePairForceTorque(int numParticles, int part1Index, int part2Index); //auxiliary

    public:
        boundary *domainBoundary;
        externalPotential<> *externalPot;
        externalPotential<vec3<double>> *externalRodPot;
        externalPotential<quaternion<double>> *externalRigidBodyPot;
        pairPotential<> *pairPot;
        pairPotential<vec3<double>, vec3<double>> *pairRodPot;
        pairPotential<quaternion<double>, quaternion<double>> *pairRigidBodyPot;
        double clock;

        /**
         * Note all potentials default to zero and not every integrator will make use of all this potentials
         * @param domainBoundary determines the boundary of the domain to be used.
         * @param externalPotential<> external potential without orientation
         * @param externalPotential<vec3<double>> external potential with rod-like orientation
         * @param pairPotential<vec3<double>,vec3<double>> pair potential between two particles with rod-like orientation
         * @param clock keeps track of global time
         */

        integrator(double dt, long seed, std::string bodytype, bool rotation);

        // Main public functions definitions
        void integrate(std::vector<particle> &parts);

        double getClock() const { return clock; }

        void setKbT(double kbt) { KbTemp = kbt; }

        // Potential and boundary related functions
        void setBoundary(boundary *bndry);

        void setExternalPotential(externalPotential<> *pot);

        void setExternalRodPotential(externalPotential<vec3<double>> *pot);

        void setExternalRigidBodyPotential(externalPotential<quaternion<double>> *pot);

        void setPairPotential(pairPotential<> *pot);

        void setPairRodPotential(pairPotential<vec3<double>, vec3<double>> *pot);

        void setPairRigidBodyPotential(pairPotential<quaternion<double>, quaternion<double>> *pot);

    };


    /**
     * Implements template functions to calculate force and torque fields
     * @tparam PARTICLE it can take values of particle or particleMS or other custom defined particles
     * @param parts is a list of PARTICLE objects
     */

    /* Calculates the total force and torque fields experienced by all particles, including external and
     * interaction pairs potentials. Relies on other functions to divide the work. */
    template< typename PARTICLE >
    void integrator::calculateTotalForceTorqueFields(std::vector<PARTICLE> &parts) {
        unsigned int N = static_cast<int>(parts.size());
        std::array<vec3<double>, 2> forctorq;
        unsigned int interactionIndex;

        // Resize force/torque fields array if the number of particles have changed
        if (static_cast<int>(forceField.size()) != N) {
            forceField.resize(N);
            torqueField.resize(N);
        }

        /* Add forces and torques due to external potential from all particles to fields,
         * initialize to zero if there is no external potential. */
        if (externalPotentialActive) {
            if (bodytype == "point") {
                for (int i = 0; i < N; i++) {
                    forctorq = externalPot->forceTorque(parts[i].position);
                    forceField[i] = forctorq[0];
                    torqueField[i] = forctorq[1];
                }
            } else if (bodytype == "rod") {
                for (int i = 0; i < N; i++) {
                    forctorq = externalRodPot->forceTorque(parts[i].position, parts[i].orientvector);
                    forceField[i] = forctorq[0];
                    torqueField[i] = forctorq[1];
                }
            } else if (bodytype == "rigidbody") {
                for (int i = 0; i < N; i++) {
                    forctorq = externalRigidBodyPot->forceTorque(parts[i].position, parts[i].orientation);
                    forceField[i] = forctorq[0];
                    torqueField[i] = forctorq[1];
                }
            } else {
                throw std::runtime_error("Unknown particle bodytype; it should be either point, rod or rigidbody.");
            }
        } else {
            for (int i = 0; i < N; i++) {
                forceField[i] = vec3<double> (0,0,0);
                torqueField[i] = vec3<double> (0,0,0);
            }
        }

        /* Add forces and torques coming from pair interactions. The routine "calculatePairsForceTorqueFields"
         * calculates all the pair forces (only once) and stores it into forcePairField and torquePairField arrays.
         * The rest of this code sums the correponding forces to the force/torque fields. Note as the
         * pair force/torque acting on i due to its interaction with j is the same as the force/torque acting on j due
         * to i (except for a minus sign), we multiply by -1 for the inverse pair interaction. The index
         * for the i,j interaction with i<j is i*(2*N - 1 - i)/2.0 + j - (i+1) (first term obtained from substracting
         * summation formulas of first N-1 and N-1-i natural numbers), with N the total number of particles.*/
        if (pairPotentialActive) {
            calculatePairsForceTorqueFields<PARTICLE>(parts);
            for (int i = 0; i < N; i++) {
                for (int j = i+1; j < N; j++) {
                    interactionIndex = i*(2*N - 1 - i)/2 + j - (i+1);
                    forceField[i] += forcePairField[interactionIndex];
                    torqueField[i] += torquePairField[interactionIndex];
                    forceField[j] += -1.0*forcePairField[interactionIndex];
                    torqueField[j] += -1.0*torquePairField[interactionIndex];
                }
            }
        }
    };


    /* Calculates all the possible pair interaction forces and torques and saves them into one dimensional
     * array of vectors: forcePairField and torquePairField. The length of this arrays corresponds to the total
     * number of possible pair interactions. The index for the i,j interaction with i<j is
     * i*(2*N - 1 - i)/2.0 + j - (i+1) (first term obtained from substracting summation formulas of first N-1 and
     * N-1-i natural numbers) with N the total number of particles. */
    template< typename PARTICLE >
    void integrator::calculatePairsForceTorqueFields(std::vector<PARTICLE> &parts) {
        unsigned int N = static_cast<int>(parts.size());
        unsigned int numberPairInteractions = N * (N - 1) / 2;
        unsigned int interactionIndex;
        std::array<vec3<double>, 2> forctorq;

        // Resize pair force/torque fields array if the number of particles have changed
        if (static_cast<int>(forcePairField.size()) != numberPairInteractions) {
            forcePairField.resize(numberPairInteractions);
            torquePairField.resize(numberPairInteractions);
        }

        // Calculate the forces and torque for each possible interaction
        if (bodytype == "point") {
            for (int i = 0; i < N; i++) {
                for (int j = i + 1; j < N; j++) {
                    forctorq = pairPot->forceTorque(parts[i].position, parts[j].position);
                    interactionIndex = i*(2*N - 1 - i)/2 + j - (i+1);
                    forcePairField[interactionIndex] = forctorq[0];
                    torquePairField[interactionIndex] = forctorq[1];
                }
            }
        } else if (bodytype == "rod") {
            for (int i = 0; i < N; i++) {
                for (int j = i + 1; j < N; j++) {
                    forctorq = pairRodPot->forceTorque(parts[i].position, parts[j].position,
                                                       parts[i].orientvector, parts[j].orientvector);
                    interactionIndex = i*(2*N - 1 - i)/2 + j - (i+1);
                    forcePairField[interactionIndex] = forctorq[0];
                    torquePairField[interactionIndex] = forctorq[1];
                }
            }
        } else if (bodytype == "rigidbody") {
            for (int i = 0; i < N; i++) {
                for (int j = i + 1; j < N; j++) {

                    forctorq = pairRigidBodyPot->forceTorque(parts[i].position, parts[j].position,
                                                             parts[i].orientation, parts[j].orientation);
                    interactionIndex = i*(2*N - 1 - i)/2 + j - (i+1);
                    forcePairField[interactionIndex] = forctorq[0];
                    torquePairField[interactionIndex] = forctorq[1];
                }
            }
        } else {
            throw std::runtime_error(
                    "Unknown particle bodytype. it should be either point, rod or rigidbody.");
        }
    };

}




