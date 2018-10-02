//
// Created by maojrs on 8/14/18.
//
#pragma once
#include "vec3.hpp"
#include "quaternion.hpp"
#include "randomgen.hpp"

namespace msmrd {
    /**
     * Abstract base class declaration for external potentials
     * @param ...AUXVARIABLES variadic template to define orientation and/or particle type, since different particle
     * types might behave different on under the same external potential. If there are no rotations or particle types
     * dependence, ...AUXVARIABLES has zero arguments); otherwise, orientation can be described by an orientation
     * vector or even a quaternion, while particle types are simply integers that provide the particle
     * type to the potential function. For example, if ...AUXVARIABLES == vec3<double> it does not depend on the
     * particle type, and its orientation describes rod-like particles, while ...AUXVARIABLES == quaternion<double>)
     * would describe general rigidbody particles with arbitrary rotations (again with no particle type dependence).
     */
    template<typename ...AUXVARIABLES>
    class externalPotential {
    public:
        externalPotential() = default;

        /* Virtual functions to calculate value of potential and force/torque at position "pos".
         * Possible orientation dependence can be added into the aux variables. */
        virtual double evaluate(vec3<double> pos, AUXVARIABLES... aux) = 0;

        virtual std::array<vec3<double>, 2> forceTorque(vec3<double> pos, AUXVARIABLES... aux) = 0;


        /* PyBind evaluation functions for external potentials of particles with no orientation, rod-like orientation
         * or full orientation. They rely on evaluate and forceTorque functions above. Note each possible combination
         * needs to be declared and implemented (template approaches failed, but might be possible). */
        double evaluatePyBind(std::vector<double> pos);
        double evaluatePyBind(std::vector<double> pos, int type);
        double evaluatePyBind(std::vector<double> pos, std::vector<double> theta);
        double evaluatePyBind(std::vector<double> pos, std::vector<double> theta, int type);

        std::vector<std::vector<double>> forceTorquePyBind(std::vector<double> pos);
        std::vector<std::vector<double>> forceTorquePyBind(std::vector<double> pos, int type);
        std::vector<std::vector<double>> forceTorquePyBind(std::vector<double> pos, std::vector<double> theta);
        std::vector<std::vector<double>> forceTorquePyBind(std::vector<double> pos, std::vector<double> theta, int type);

    };


    /**
     * Abstract base class declaration for pair potentials
     * @param ...AUXVARIABLES variadic template to define orientation and/or particle type. If there are no
     * rotations nor particle types, ...AUXVARIABLES has zero arguments); otherwise, orientation can be described by
     * an orientation vector or even a quaternion, while particle types are simply integers that provide the particle
     * type to the potential function. For example, if ...AUXVARIABLES == vec3<double>, vec3<double> it does not
     * depend on the particle type, and its orientation describes pairs of rod-like particles,
     * while ...AUXVARIABLES == quaternion<double>, quaternion<double> would describe general rigidbody particles
     * with arbitrary rotations (again with no particle type dependence).
     */
    template<typename ...AUXVARIABLES>
    class pairPotential {
    public:
        pairPotential() = default;

        /* Virtual functions to calculate value of potential and force/torque at positions "pos1" and "pos2".
         * Possible orientation dependence can be added into the aux variables. The function forceTorque should
         * return (force1, torque1, force2, torque2), the first two correspond to the force and torque acting on
         * particle 1 due to its interaction with particle 2, and the second two correspond to the force and
         * torque acting on particle 2 due to its interaction with particle 1.*/
        virtual double evaluate(vec3<double> pos1, vec3<double> pos2, AUXVARIABLES... aux) = 0;

        virtual std::array<vec3<double>, 4> forceTorque(vec3<double> pos1, vec3<double> pos2, AUXVARIABLES... aux) = 0;


        /* PyBind evaluation functions for pair potentials of particles with no orientation, rod-like orientation or
         * full orientation. They rely on evaluate and forceTorque functions above. Note each possible combination
         * needs to be declared and implemented (template approaches failed, but might be possible). */
        double evaluatePyBind(std::vector<double> pos1, std::vector<double> pos2);
        double evaluatePyBind(std::vector<double> pos1, std::vector<double> pos2, int type1, int type2);
        double evaluatePyBind(std::vector<double> pos1, std::vector<double> pos2,
                              std::vector<double> theta1, std::vector<double> theta2);
        double evaluatePyBind(std::vector<double> pos1, std::vector<double> pos2,
                              std::vector<double> theta1, std::vector<double> theta2,
                              int type1, int type2);

        std::vector<std::vector<double>> forceTorquePyBind(std::vector<double> pos1, std::vector<double> pos2);
        std::vector<std::vector<double>> forceTorquePyBind(std::vector<double> pos1, std::vector<double> pos2,
                                                           int type1, int type2);
        std::vector<std::vector<double>> forceTorquePyBind(std::vector<double> pos1, std::vector<double> pos2,
                                                           std::vector<double> theta1, std::vector<double> theta2);
        std::vector<std::vector<double>> forceTorquePyBind(std::vector<double> pos1, std::vector<double> pos2,
                                                           std::vector<double> theta1, std::vector<double> theta2,
                                                           int type1, int type2);

    };


    template<long unsigned int N>
    std::vector<std::vector<double>> array2vec(std::array<vec3<double>, N> forceTorque) {
        std::vector<std::vector<double>> output;
        output.resize(N);
        for (int i = 0; i < N; i++) {
            output[i].resize(3);
            output[i][0] = 1.0*forceTorque[i][0];
            output[i][1] = 1.0*forceTorque[i][1];
            output[i][2] = 1.0*forceTorque[i][2];
        }
        return output;
    }
}


