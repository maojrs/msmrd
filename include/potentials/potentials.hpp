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
     * @param ...ORIENTATION variadic template to define orientation. If there are no rotations, ...ORIENTATION has
     * zero arguments); otherwise, it can be described by an orientation vector or even a quaternion. If the
     * orientation is simply a vector (vec3<double>) it describes rod-like particles.
     */
    template<typename ...ORIENTATION>
    class externalPotential {
    public:
        externalPotential() = default;

        // Virtual functions to calculate value of potential and force/torque at position "pos" and orientation u
        virtual double evaluate(vec3<double> pos, ORIENTATION... theta) = 0;
        virtual std::array<vec3<double>, 2> forceTorque(vec3<double> pos, ORIENTATION... theta) = 0;


        /* PyBind evaluation functions for external potentials of particles with no orientation, rod-like orientation or
         * full orientation. They rely on evaluate and forceTorque functions above */
        double evaluatePyBind(std::vector<double> pos);
        double evaluatePyBind(std::vector<double> pos, std::vector<double> theta);
        std::vector<double> forceTorquePyBind(std::vector<double> pos);
        std::vector<std::vector<double>> forceTorquePyBind(std::vector<double> pos, std::vector<double> theta);

    };


    /**
     * Abstract base class declaration for pair potentials
     * @param ...ORIENTATION variadic template to define orientation. If there are no rotations, ...ORIENTATION has
     * zero arguments); otherwise, it can be described by an orientation vector or even a quaternion. If the
     * orientation is simply a vector (vec3<double>) it describes rod-like particles.
     */
    template<typename ...ORIENTATION>
    class pairPotential {
    public:
        pairPotential() = default;

        /* Virtual functions to calculate value of potential and force/torque at positions "pos1" and "pos2"
         * and orientations u. The function forceTorque should return (force1, torque1, force2, torque2), which
         * correspond to the force and torque acting on particle 1 and 2 respectively.*/
        virtual double evaluate(vec3<double> pos1, vec3<double> pos2, ORIENTATION... theta) = 0;
        virtual std::array<vec3<double>, 4> forceTorque(vec3<double> pos1, vec3<double> pos2, ORIENTATION... theta) = 0;


        /* PyBind evaluation functions for pair potentials of particles with no orientation, rod-like orientation or
         * full orientation. They rely on evaluate and forceTorque functions above */
        double evaluatePyBind(std::vector<double> pos1, std::vector<double> pos2);
        double evaluatePyBind(std::vector<double> pos1, std::vector<double> pos2,
                              std::vector<double> theta1, std::vector<double> theta2);
        std::vector<double> forceTorquePyBind(std::vector<double> pos1, std::vector<double> pos2);
        std::vector<std::vector<double>> forceTorquePyBind(std::vector<double> pos1,
                                                           std::vector<double> pos2,
                                                           std::vector<double> theta1,
                                                           std::vector<double> theta2);



    };

}


