//
// Created by maojrs on 8/14/18.
//
#pragma once
#include "particle.hpp"
#include "randomgen.hpp"
#include "quaternion.hpp"
#include "vec3.hpp"

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
        virtual double evaluate(const particle &part) = 0;

        virtual std::array<vec3<double>, 2> forceTorque(const particle &part) = 0;

        std::vector<std::vector<double>> forceTorquePyBind(const particle &part);
    };


    /**
     * Abstract base class declaration for pair potentials (note particle type can take also child of
     * particle class like particleMS)
     */
    class pairPotential {
    public:
        pairPotential() = default;

        /* Virtual functions to calculate value of potential and force/torque at positions "pos1" and "pos2".
         * Possible orientation dependence can be added into the aux variables. The function forceTorque should
         * return (force1, torque1, force2, torque2), the first two correspond to the force and torque acting on
         * particle 1 due to its interaction with particle 2, and the second two correspond to the force and
         * torque acting on particle 2 due to its interaction with particle 1.*/
        virtual double evaluate(const particle &part1, const particle &part2) = 0;

        virtual std::array<vec3<double>, 4> forceTorque(const particle &part1, const particle &part2) = 0;

        std::vector<std::vector<double>> forceTorquePyBind(const particle &part1, const particle &part2);

    };


}


