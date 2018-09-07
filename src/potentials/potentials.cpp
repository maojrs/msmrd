//
// Created by maojrs on 8/21/18.
//
#include "vec3.hpp"
#include "potentials/potentials.hpp"

namespace msmrd {
    /**
     *  Implementation of non-abstract functions of externalPotential abstract class (default constructor in header)
     */

    /* Needed for PyBinding externalPotential.evaluate function , since it can only take vectors as input.
     * Note the template ORIENTATION takes zero arguments, so it corresponds to no orientation */
    template<>
    double externalPotential<>::evaluatePyBind(std::vector<double> pos) {
        vec3<double> x = vec3<double>(pos);
        return evaluate(x);
    }

    /* Needed for PyBinding externalPotential.forceTorque function , since it can only take vectors as input.
     * Note the template ORIENTATION takes zero arguments, so it corresponds to no orientation */
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

    /* Needed for PyBinding externalPotential.evaluate function of rod-like particles. Note the template
     * ORIENTATION corresponds to vec3<double>, so it corresponds to rod-like particles with orientation
     * described by one vector. */
    template<>
    double externalPotential<vec3<double>>::evaluatePyBind(std::vector<double> pos1, std::vector<double> u) {
        vec3<double> x = vec3<double>(pos1);
        vec3<double> theta = vec3<double>(u);
        return evaluate(x, theta);
    }

    /* Needed for PyBinding externalPotential.forceTorque function of rod-like particles. Note the template
     * ORIENTATION corresponds to vec3<double>, so it corresponds to rod-like particles with orientation
     * described by one vector. */
    template<>
    std::vector<std::vector<double>>
    externalPotential<vec3<double>>::forceTorquePyBind(std::vector<double> pos1, std::vector<double> u) {
        vec3<double> x = vec3<double>(pos1);
        vec3<double> theta = vec3<double>(u);
        std::array<vec3<double>, 2> forceTorquex = forceTorque(x, theta);
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
     */


    /* Needed for PyBinding pairPotential.evaluate function. Note the template ORIENTATION takes
     * zero arguments, so it corresponds to no orientation. */
    template<>
    double pairPotential<>::evaluatePyBind(std::vector<double> pos1, std::vector<double> pos2) {
        vec3<double> x1 = vec3<double>(pos1);
        vec3<double> x2 = vec3<double>(pos2);
        return evaluate(x1, x2);
    }

    /* Needed for PyBinding pairPotential.forceTorque function , since it can only take vectors as input.
     * Note the template ORIENTATION takes zero arguments, so it corresponds to no orientation */
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


    /* Needed for PyBinding pairPotential.evaluate function of rod-like particles. Note the template
     * ORIENTATION corresponds to <vec3<double>,vec3<double>>, so it corresponds to two rod-like particles each
     * with orientation described by one vector. */
    template<>
    double pairPotential<vec3<double>, vec3<double>>::evaluatePyBind(std::vector<double> pos1,
                                                                     std::vector<double> pos2,
                                                                     std::vector<double> u1,
                                                                     std::vector<double> u2) {
        vec3<double> x1 = vec3<double>(pos1);
        vec3<double> x2 = vec3<double>(pos2);
        vec3<double> theta1 = vec3<double>(u1);
        vec3<double> theta2 = vec3<double>(u2);
        return evaluate(x1, x2, theta1, theta2);
    }

    /* Needed for PyBinding pairPotential.forceTorque function of rod-like particles. Note the template
     * ORIENTATION corresponds to <vec3<double>,vec3<double>>, so it corresponds to two rod-like particles each
     * with orientation described by one vector. */
    template<>
    std::vector<std::vector<double>>
    pairPotential<vec3<double>, vec3<double>>::forceTorquePyBind(std::vector<double> pos1,
                                                                 std::vector<double> pos2,
                                                                 std::vector<double> u1,
                                                                 std::vector<double> u2) {
        vec3<double> x1 = vec3<double>(pos1);
        vec3<double> x2 = vec3<double>(pos2);
        vec3<double> theta1 = vec3<double>(u1);
        vec3<double> theta2 = vec3<double>(u2);
        std::array<vec3<double>, 2> forceTorquex = forceTorque(x1, x2, theta1, theta2);
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