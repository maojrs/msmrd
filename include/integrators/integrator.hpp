//
// Created by dibakma on 27.06.18.
//

#pragma once
#include <array>
#include <utility>
#include <memory>
#include "boundary.hpp"
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
        randomgen randg;
        bool rotation;
        bool boundaryActive = false;
        bool externalPotentialActive = false;
        bool pairPotentialActive = false;


        /**
        * @param KbTemp = Boltzman constant times temperature
        * @param dt time step
        * @param seed variable for random number generation (Note seed = -1 corresponds to random device)
        * @param randg random number generator based in mt19937
        * @param rotation boolean to indicate if rotation should be integrated
        * @param externalPotActive indicates if external potential has been set
        * @param pairPotActive indicates if potential potential has been set
        */

        // Protected abstract functions
        virtual void integrateOne(int index, std::vector<particle> &parts,
                                  double timestep) = 0; // Version for pair interactions (full particle list required)
        virtual void translate(particle &part, vec3<double> force, double dt) = 0;

        virtual void rotate(particle &part, vec3<double> torque, double dt) = 0;

        // Protected functions to get forces and torques due to external or pair potentials for integrator
        std::array<vec3<double>, 2> getExternalForceTorque(particle &part);

        std::array<vec3<double>, 2> getPairsForceTorque(int partIndex, std::vector<particle> &parts);

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

        integrator(double dt, long seed, bool rotation);

        // Main functions definitions
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

//    double evalExternalPotential(std::vector<double> pos);
//    double evalPairPotential(std::vector<double> pos1, std::vector<double> pos2);
//    double evalRodPairPotential(std::vector<double> pos1,
//                                std::vector<double> pos2,
//                                std::vector<double> u1,
//                                std::vector<double> u2);
//    std::vector<double> evalExternalForce(std::vector<double> pos);
//    std::vector<double> evalPairForce(std::vector<double> pos1, std::vector<double> pos2);
//    std::vector<std::vector<double>> evalRodPairForce(std::vector<double> pos1,
//                                                      std::vector<double> pos2,
//                                                      std::vector<double> u1,
//                                                      std::vector<double> u2);
    };

}







