//
// Created by maojrs on 8/21/18.
//
#include "vec3.hpp"
#include "potentials/potentials.hpp"

namespace msmrd {
    /**
     *  Implementation of non-abstract functions of externalPotential abstract class (default constructor in header)
     *  All these functions are needed for PyBinding evaluate and forceTorque functions, since pyBind only accepts
     *  vectors as input. The different versions of the same functions depend on the ORIENTATION template that can
     *  be used for particles with no orientation, rod-like particles with their orientation described by one vectors,
     *  and by rigidbody particles with no symmetry and their orientation described by quaternions.
     */


    //Evaluates potential. Template ORIENTATION takes zero arguments, so it corresponds to no orientation.
    template<>
    double externalPotential<>::evaluatePyBind(std::vector<double> pos) {
        vec3<double> x = vec3<double>(pos);
        return evaluate(x);
    }


    /* Evaluates potential. Template ORIENTATION corresponds to vec3<double>, so it corresponds to rod-like
     * particles with orientation described by one vector. */
    template<>
    double externalPotential<vec3<double>>::evaluatePyBind(std::vector<double> pos1, std::vector<double> theta1) {
        vec3<double> x = vec3<double>(pos1);
        vec3<double> th = vec3<double>(theta1);
        return evaluate(x, th);
    }


    /* Evaluates potential. Template ORIENTATION corresponds to quaternion<double>, so it corresponds to
     * rigidbody particles with orientation described by a quaternion. */
    template<>
    double externalPotential<quaternion<double>>::evaluatePyBind(std::vector<double> pos1, std::vector<double> theta1) {
        vec3<double> x = vec3<double>(pos1);
        quaternion<double> th = quaternion<double>(theta1);
        return evaluate(x, th);
    }


    // Evaluates force and torque. Template ORIENTATION takes zero arguments, so it corresponds to no orientation.
    template<>
    std::vector<double> externalPotential<>::forceTorquePyBind(std::vector<double> pos) {
        vec3<double> x = vec3<double>(pos);
        std::array<vec3<double>, 2> forcex = forceTorque(x);
        std::vector<double> output;
        output.resize(3);
        output[0] = forcex[0][0];
        output[1] = forcex[0][1];
        output[2] = forcex[0][2];
        return output;
    }


    /* Evaluates force and torque. Template ORIENTATION corresponds to vec3<double>, so it corresponds to
     * rod-like particles with orientation described by one vector. */
    template<>
    std::vector<std::vector<double>>
    externalPotential<vec3<double>>::forceTorquePyBind(std::vector<double> pos1, std::vector<double> theta1) {
        vec3<double> x = vec3<double>(pos1);
        vec3<double> th = vec3<double>(theta1);
        std::array<vec3<double>, 2> forceTorquex = forceTorque(x, th);
        std::vector<std::vector<double>> output;
        output.resize(2);
        output[0].resize(3);
        output[1].resize(3);
        for (int i = 0; i < 2; i++) {
            output[i][0] = forceTorquex[i][0];
            output[i][1] = forceTorquex[i][1];
            output[i][2] = forceTorquex[i][2];
        }
        return output;
    }


    /* Evaluates force and torque. Template ORIENTATION corresponds to quaternion<double>, so it corresponds to
     * rigidbody particles with orientation described by a quaternion. */
    template<>
    std::vector<std::vector<double>>
    externalPotential<quaternion<double>>::forceTorquePyBind(std::vector<double> pos1, std::vector<double> theta1) {
        vec3<double> x = vec3<double>(pos1);
        quaternion<double> th = quaternion<double>(theta1);
        std::array<vec3<double>, 2> forceTorquex = forceTorque(x, th);
        std::vector<std::vector<double>> output;
        output.resize(2);
        output[0].resize(3);
        output[1].resize(3);
        for (int i = 0; i < 2; i++) {
            output[i][0] = forceTorquex[i][0];
            output[i][1] = forceTorquex[i][1];
            output[i][2] = forceTorquex[i][2];
        }
        return output;
    }


    /**
     *  Implementation of non-abstract functions of pairPotential abstract class (default constructor in header)
     *  All these functions are needed for PyBinding evaluate and forceTorque functions, since pyBind only accepts
     *  vectors as input. The different versions of the same functions depend on the ORIENTATION template that can
     *  be used for particles with no orientation, rod-like particles with their orientation described by one vectors,
     *  and by rigidbody particles with no symmetry and their orientation described by quaternions.
     */


    // Evaluates potential. Template ORIENTATION takes zero arguments, so it corresponds to no orientation.
    template<>
    double pairPotential<>::evaluatePyBind(std::vector<double> pos1, std::vector<double> pos2) {
        vec3<double> x1 = vec3<double>(pos1);
        vec3<double> x2 = vec3<double>(pos2);
        return evaluate(x1, x2);
    }


    /* Evaluate potential. Template ORIENTATION corresponds to <vec3<double>, vec3<double>>, so it
     * corresponds to two rod-like particles each with orientation described by one vector. */
    template<>
    double pairPotential<vec3<double>, vec3<double>>::evaluatePyBind(std::vector<double> pos1,
                                                                     std::vector<double> pos2,
                                                                     std::vector<double> theta1,
                                                                     std::vector<double> theta2) {
        vec3<double> x1 = vec3<double>(pos1);
        vec3<double> x2 = vec3<double>(pos2);
        vec3<double> th1 = vec3<double>(theta1);
        vec3<double> th2 = vec3<double>(theta2);
        return evaluate(x1, x2, th1, th2);
    }


    /* Evaluate potential. Template ORIENTATION corresponds to <quaternion<double>, quaternion<double>>, so it
     * corresponds to two rigidbody particles each with orientation described by a quaternion. */
    template<>
    double pairPotential<quaternion<double>, quaternion<double>>::evaluatePyBind(std::vector<double> pos1,
                                                                                 std::vector<double> pos2,
                                                                                 std::vector<double> theta1,
                                                                                 std::vector<double> theta2) {
        vec3<double> x1 = vec3<double>(pos1);
        vec3<double> x2 = vec3<double>(pos2);
        quaternion<double> th1 = quaternion<double>(theta1);
        quaternion<double> th2 = quaternion<double>(theta2);
        return evaluate(x1, x2, th1, th2);
    }


    // Evaluates force and torque. Template ORIENTATION takes zero arguments, so it corresponds to no orientation
    template<>
    std::vector<double> pairPotential<>::forceTorquePyBind(std::vector<double> pos1, std::vector<double> pos2) {
        vec3<double> x1 = vec3<double>(pos1);
        vec3<double> x2 = vec3<double>(pos2);
        std::array<vec3<double>, 2> forcex = forceTorque(x1, x2);
        std::vector<double> output;
        output.resize(3);
        output[0] = forcex[0][0];
        output[1] = forcex[0][1];
        output[2] = forcex[0][2];
        return output;
    }


    /* Evaluates force and torque. Template ORIENTATION corresponds to <vec3<double>,vec3<double>>, so
     * it corresponds to two rod-like particles each with orientation described by one vector. */
    template<>
    std::vector<std::vector<double>>
    pairPotential<vec3<double>, vec3<double>>::forceTorquePyBind(std::vector<double> pos1,
                                                                 std::vector<double> pos2,
                                                                 std::vector<double> theta1,
                                                                 std::vector<double> theta2) {
        vec3<double> x1 = vec3<double>(pos1);
        vec3<double> x2 = vec3<double>(pos2);
        vec3<double> th1 = vec3<double>(theta1);
        vec3<double> th2 = vec3<double>(theta2);
        std::array<vec3<double>, 2> forceTorquex = forceTorque(x1, x2, th1, th2);
        std::vector<std::vector<double>> output;
        output.resize(2);
        output[0].resize(3);
        output[1].resize(3);
        for (int i = 0; i < 2; i++) {
            output[i][0] = forceTorquex[i][0];
            output[i][1] = forceTorquex[i][1];
            output[i][2] = forceTorquex[i][2];
        }
        return output;
    }


    /* Evaluates force and torque. Template ORIENTATION corresponds to <quaternion<double>, quaternion<double>>,
     * so it corresponds to two rigidbody particles each with orientation described by a quaternion. */
    template<>
    std::vector<std::vector<double>>
    pairPotential<quaternion<double>, quaternion<double>>::forceTorquePyBind(std::vector<double> pos1,
                                                                             std::vector<double> pos2,
                                                                             std::vector<double> theta1,
                                                                             std::vector<double> theta2) {
        vec3<double> x1 = vec3<double>(pos1);
        vec3<double> x2 = vec3<double>(pos2);
        quaternion<double> th1 = quaternion<double>(theta1);
        quaternion<double> th2 = quaternion<double>(theta2);
        std::array<vec3<double>, 2> forceTorquex = forceTorque(x1, x2, th1, th2);
        std::vector<std::vector<double>> output;
        output.resize(2);
        output[0].resize(3);
        output[1].resize(3);
        for (int i = 0; i < 2; i++) {
            output[i][0] = forceTorquex[i][0];
            output[i][1] = forceTorquex[i][1];
            output[i][2] = forceTorquex[i][2];
        }
        return output;
    }

}