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
         * or full orientation. They rely on evaluate and forceTorque functions above */

//        template<typename ...AUXVARS>
//        double evaluatePyBind(std::vector<double> pos, AUXVARS... auxvars);
//        template<typename ...AUXVARS>
//        std::vector<std::vector<double>> forceTorquePyBind(std::vector<double> pos, AUXVARS... auxvars);

        double evaluatePyBind(std::vector<double> pos);
        double evaluatePyBind(std::vector<double> pos, std::vector<double> theta);
        double evaluatePyBind(std::vector<double> pos, std::vector<double> theta, int type);

        std::vector<double> forceTorquePyBind(std::vector<double> pos);
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
         * full orientation. They rely on evaluate and forceTorque functions above */
//        template<typename ...AUXVARS>
//        double evaluatePyBind(std::vector<double> pos1, std::vector<double> pos2, AUXVARS... auxvars);
//        template<typename ...AUXVARS>
//        std::vector<std::vector<double>> forceTorquePyBind(std::vector<double> pos1, std::vector<double> pos2,
//                                                           AUXVARS... auxvars);


        double evaluatePyBind(std::vector<double> pos1, std::vector<double> pos2);
        double evaluatePyBind(std::vector<double> pos1, std::vector<double> pos2,
                              std::vector<double> theta1, std::vector<double> theta2);

        std::vector<double> forceTorquePyBind(std::vector<double> pos1, std::vector<double> pos2);
        std::vector<std::vector<double>> forceTorquePyBind(std::vector<double> pos1, std::vector<double> pos2,
                                                           std::vector<double> theta1, std::vector<double> theta2);

    };





//    //Evaluates potential. Template AUXVARIABLES takes zero arguments, so it corresponds to no orientation.
//    template<typename ...AUXVARIABLES>
//    template<typename ...AUXVARS>
//    double externalPotential<AUXVARIABLES...>::evaluatePyBind(std::vector<double> pos1, AUXVARS... auxvars) {
//        int n = sizeof...(AUXVARS);
//        std::tuple<AUXVARS...> auxvarsPack(auxvars...); // Stores auxvars into auxvarsPack tuple
//
//        vec3<double> x = vec3<double>(pos1);
//
//        if ( n == 0 ) {
//            return evaluate(x);
//        }
//        else if (n >= 1 && n < 3) {
//            // extract type of first argument in variadic template
//            using ORIENTATION = std::tuple_element<0, std::tuple<AUXVARIABLES...>>;
//            ORIENTATION theta1 = ORIENTATION(auxvarsPack[0]);
//            if ( n == 1 ) { return evaluate(x, theta1); };
//            if ( n == 2 ) { return evaluate(x, theta1, std::get<1>(auxvarsPack)); };
//        }
//        else{
//            throw std::runtime_error("Variadic template for external potential only supports up to two arguments.");
//        }
//    }
//
//    // Evaluates force and torque. Template AUXVARIABLES takes zero arguments, so it corresponds to no orientation.
//    template<typename ...AUXVARIABLES>
//    template<typename ...AUXVARS>
//    std::vector<std::vector<double>> externalPotential<AUXVARIABLES...>::forceTorquePyBind(std::vector<double> pos1,
//                                                                                           AUXVARS... auxvars) {
//        int n = sizeof...(AUXVARS);
//        std::tuple<AUXVARS...> auxvarsPack(auxvars...); // Stores auxvars into auxvarsPack tuple
//
//        vec3<double> x = vec3<double>(pos1);
//        std::array<vec3<double>, 2> forcetorq;
//        std::vector<std::vector<double>> output;
//
//        if ( n == 0 ) {
//            forcetorq = forceTorque(x, auxvarsPack);
//        }
//        else if (n >= 1 && n < 3) {
//            // extract type of first argument in variadic template
//            using ORIENTATION = typename std::tuple_element<0, std::tuple<AUXVARIABLES...>>;
//            ORIENTATION theta1 = ORIENTATION(auxvarsPack[0]);
//            if ( n == 1 ) { forcetorq = forceTorque(x, theta1); };
//            if ( n == 2 ) { forcetorq = forceTorque(x, theta1, auxvarsPack[1]); };
//        }
//        else{
//            throw std::runtime_error("Variadic template for external potential only supports up to two arguments.");
//        }
//
//        output.resize(2);
//        output[0].resize(3);
//        output[1].resize(3);
//        for (int i = 0; i < 2; i++) {
//            output[i][0] = 1.0*forcetorq[i][0];
//            output[i][1] = 1.0*forcetorq[i][1];
//            output[i][2] = 1.0*forcetorq[i][2];
//        }
//        return output;
//    }
//
//
//    //Evaluates potential. Template AUXVARIABLES takes zero arguments, so it corresponds to no orientation.
//    template<typename ...AUXVARIABLES>
//    template<typename ...AUXVARS>
//    double pairPotential<AUXVARIABLES...>::evaluatePyBind(std::vector<double> pos1, std::vector<double> pos2,
//                                                          AUXVARS... auxvars) {
//        int n = sizeof...(AUXVARS);
//        std::tuple<AUXVARS...> auxvarsPack(auxvars...); // Stores auxvars into auxvarsPack tuple
//
//        vec3<double> x1 = vec3<double>(pos1);
//        vec3<double> x2 = vec3<double>(pos2);
//        //std::tuple<AUXVARS...> auxvariables(auxvars...);
//        if ( n == 0 ) {
//            return evaluate(x1, x2, auxvarsPack);
//        }
//        else if (n == 2 || n == 4) {
//            // extract type of first argument in variadic template
//            using ORIENTATION = typename std::tuple_element<0, std::tuple<AUXVARIABLES...>>;
//            ORIENTATION theta1 = ORIENTATION(auxvarsPack[0]);
//            ORIENTATION theta2 = ORIENTATION(auxvarsPack[1]);
//            if ( n == 2 ) { return evaluate(x1, x2, theta1, theta2); };
//            if ( n == 4 ) { return evaluate(x1, x2, theta1, theta2, auxvarsPack[2], auxvarsPack[3]); };
//        }
//        else{
//            throw std::runtime_error("Variadic template for pair potential only supports zero, two or four arguments.");
//        }
//    }
//
//
//    // Evaluates force and torque. Template AUXVARIABLES takes zero arguments, so it corresponds to no orientation.
//    template<typename ...AUXVARIABLES>
//    template<typename ...AUXVARS>
//    std::vector<std::vector<double>> pairPotential<AUXVARIABLES...>::forceTorquePyBind(std::vector<double> pos1,
//                                                                                       std::vector<double> pos2,
//                                                                                       AUXVARS... auxvars) {
//        int n = sizeof...(AUXVARS);
//        std::tuple<AUXVARS...> auxvarsPack(auxvars...); // Stores auxvars into auxvarsPack tuple
//
//        vec3<double> x1 = vec3<double>(pos1);
//        vec3<double> x2 = vec3<double>(pos2);
//        std::array<vec3<double>, 4> forcetorq;
//        std::vector<std::vector<double>> output;
//
//        if ( n == 0 ) {
//            forcetorq = forceTorque(x1, x2, auxvarsPack);
//        }
//        else if (n == 2 || n == 4) {
//            // extract type of first argument in variadic template
//            using ORIENTATION = typename std::tuple_element<0, std::tuple<AUXVARIABLES...>>;
//            ORIENTATION theta1 = ORIENTATION(auxvarsPack[0]);
//            ORIENTATION theta2 = ORIENTATION(auxvarsPack[1]);
//            if ( n == 2 ) { forcetorq = forceTorque(x1, x2, theta1, theta2); };
//            if ( n == 4 ) { forcetorq = forceTorque(x1, x2, theta1, theta2, auxvarsPack[2], auxvarsPack[3]); };
//        }
//        else{
//            throw std::runtime_error("Variadic template for pair potential only supports zero, two or four arguments.");
//        }
//
//        output.resize(4);
//        for (int i = 0; i < 4; i++) {
//            output[i].resize(3);
//            output[i][0] = 1.0*forcetorq[i][0];
//            output[i][1] = 1.0*forcetorq[i][1];
//            output[i][2] = 1.0*forcetorq[i][2];
//        }
//        return output;
//    }

}


