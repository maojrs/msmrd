//
// Created by maojrs on 8/14/18.
//
#pragma once
#include <memory>
#include "particle.hpp"
#include "randomgen.hpp"
#include "quaternion.hpp"
#include "tools.hpp"
#include "vec3.hpp"
#include "boundaries/boundary.hpp"

namespace msmrd {
    /**
     * Abstract base class declaration for external potentials (note particle type can take also child of
     * particle class like particleMS)
     */
    class externalPotential {
    public:
        externalPotential() = default;

        /* Virtual functions to calculate value of potential and force/torque at position "pos".
         * Possible orientation dependence can be added into the aux variables. */
        virtual double evaluate(particle &part) = 0;

        virtual std::array<vec3<double>, 2> forceTorque(particle &part) = 0;

        std::vector<std::vector<double>> forceTorquePyBind(particle &part);
    };


    /**
     * Abstract base class declaration for pair potentials (note particle type can take also child of
     * particle class like particleMS)
     */
    class pairPotential {
    protected:
        boundary *domainBoundary;
        bool boundaryActive = false;
        /*
        * @param *domainBoundary pointer to the boundary object to be used. Useful to compute potentials
        * in periodic domains. It must point to the same boundary as the integrator.
        * @param boundaryActive true is boundary is active in the system.
        */
    public:

        pairPotential() = default;

        /* Virtual functions to calculate value of potential and force/torque at positions "pos1" and "pos2".
         * Possible orientation dependence can be added into the aux variables. The function forceTorque should
         * return (force1, torque1, force2, torque2), the first two correspond to the force and torque acting on
         * particle 1 due to its interaction with particle 2, and the second two correspond to the force and
         * torque acting on particle 2 due to its interaction with particle 1.*/
        virtual double evaluate(particle &part1, particle &part2) = 0;

        virtual std::array<vec3<double>, 4> forceTorque(particle &part1, particle &part2) = 0;


        // Function to translate forceTorque function to pyBind
        std::vector<std::vector<double>> forceTorquePyBind(particle &part1, particle &part2);


        // Additional useful functions so potentials can deal with possible periodic boundaries.
        void setBoundary(boundary *bndry);

        vec3<double> relativePosition(const vec3<double> p1, const vec3<double> p2);

        std::array<vec3<double>, 2> relativePositionComplete(vec3<double> p1, vec3<double> p2);

    };


    /**
     * Class to combine several pair potentials into one. Potentials need to be added with the addPotential
     * function.
     */
     class combinedPairPotential : public pairPotential{
     public:
         std::vector<std::shared_ptr<pairPotential>> potentials;
         /*
          * @param pairPotentials, tuple of pointers to pairPotentials. The evaluate and forceTorque
          * functions call all the potentials in the tuple to evaluate the functions.
          */

         using pairPotential::pairPotential;

         void addPotential(std::shared_ptr<pairPotential> pairPotentialPtr);

         double evaluate(particle &part1, particle &part2) override;

         std::array<vec3<double>, 4> forceTorque(particle &part1, particle &part2) override;
     };

    /**
    * Class to combine several external potentials into one. Potentials need to be added with the
     * addPotential function.
    */
    class combinedExternalPotential : public externalPotential{
    public:
        std::vector<std::shared_ptr<externalPotential>> potentials;
        /*
         * @param pairPotentials, tuple of pointers to pairPotentials. The evaluate and forceTorque
         * functions call all the potentials in the tuple to evaluate the functions.
         */

        using externalPotential::externalPotential;

        void addPotential(std::shared_ptr<externalPotential> externalPotentialPtr);

        double evaluate(particle &part) override;

        std::array<vec3<double>, 2> forceTorque(particle &part) override;
    };


}


