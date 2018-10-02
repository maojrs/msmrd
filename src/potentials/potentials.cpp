//
// Created by maojrs on 8/21/18.
//
#include <tuple>
#include "potentials/potentials.hpp"
#include "tools.hpp"
#include "vec3.hpp"

namespace msmrd {
    /**
     *  Implementation of non-abstract functions of externalPotential abstract class (default constructor in header)
     *  All these functions are needed for PyBinding evaluate and forceTorque functions, since pyBind only accepts
     *  vectors as input. The different versions of the same functions depend on the AUXVARIABLES template that can
     *  be used for particles with no orientation, rod-like particles with their orientation described by one vectors,
     *  and by rigidbody particles with no symmetry and their orientation described by quaternions. AUXVARIABLES
     *  can also incorporate particle types to incorporate type-specific behavior in the potential function.
     */


    /* Evaluate functions */


    /* Evaluates potential. Template AUXVARIABLES is empty, so it corresponds to particles with
     * no relevant orientation. */
    template<>
    double externalPotential<>::evaluatePyBind(std::vector<double> pos1) {
        vec3<double> x = vec3<double>(pos1);
        return evaluate(x);
    }

    /* Same as previous function with particle type dependence */
    template<>
    double externalPotential<int>::evaluatePyBind(std::vector<double> pos1, int type1) {
        vec3<double> x = vec3<double>(pos1);
        return evaluate(x, type1);
    }

    /* Evaluates potential. Template AUXVARIABLES corresponds to vec3<double>, so it corresponds to rod-like
     * particles with orientation described by one vector. */
    template<>
    double externalPotential<vec3<double>>::evaluatePyBind(std::vector<double> pos1, std::vector<double> theta1) {
        vec3<double> x = vec3<double>(pos1);
        vec3<double> th = vec3<double>(theta1);
        return evaluate(x, th);
    }

    /* Same as previous function with particle type dependence */
    template<>
    double externalPotential<vec3<double>, int>::evaluatePyBind(std::vector<double> pos1, std::vector<double> theta1,
                                                                int type1) {
        vec3<double> x = vec3<double>(pos1);
        vec3<double> th = vec3<double>(theta1);
        return evaluate(x, th, type1);
    }


    /* Evaluates potential. Template AUXVARIABLES corresponds to quaternion<double>, so it corresponds to
     * rigidbody particles with orientation described by a quaternion. */
    template<>
    double externalPotential<quaternion<double>>::evaluatePyBind(std::vector<double> pos1, std::vector<double> theta1) {
        vec3<double> x = vec3<double>(pos1);
        quaternion<double> th = quaternion<double>(theta1);
        return evaluate(x, th);
    }

    /* Same as previous function with particle type dependence */
    template<>
    double externalPotential<quaternion<double>, int>::evaluatePyBind(std::vector<double> pos1,
                                                                      std::vector<double> theta1, int type1) {
        vec3<double> x = vec3<double>(pos1);
        quaternion<double> th = quaternion<double>(theta1);
        return evaluate(x, th, type1);
    }


    /* Force torque functions */

    /* Evaluates force and torque. Template AUXVARIABLES takes zero arguments, so it corresponds to no orientation. */
    template<>
    std::vector<std::vector<double>> externalPotential<>::forceTorquePyBind(std::vector<double> pos) {
        vec3<double> x = vec3<double>(pos);
        std::array<vec3<double>, 2> forceTorquex = forceTorque(x);
        return msmrdtools::array2Dtovec2D(forceTorquex);
    }

    /* Same as previous function with particle type dependence */
    template<>
    std::vector<std::vector<double>> externalPotential<int>::forceTorquePyBind(std::vector<double> pos1, int type1) {
        vec3<double> x = vec3<double>(pos1);
        std::array<vec3<double>, 2> forceTorquex = forceTorque(x, type1);
        return msmrdtools::array2Dtovec2D(forceTorquex);
    }


    /* Evaluates force and torque. Template AUXVARIABLES corresponds to vec3<double>, so it corresponds to
     * rod-like particles with orientation described by one vector. */
    template<>
    std::vector<std::vector<double>>
    externalPotential<vec3<double>>::forceTorquePyBind(std::vector<double> pos1, std::vector<double> theta1) {
        vec3<double> x = vec3<double>(pos1);
        vec3<double> th = vec3<double>(theta1);
        std::array<vec3<double>, 2> forceTorquex = forceTorque(x, th);
        return msmrdtools::array2Dtovec2D(forceTorquex);
    }

    /* Same as previous function with particle type dependence */
    template<>
    std::vector<std::vector<double>>
    externalPotential<vec3<double>, int>::forceTorquePyBind(std::vector<double> pos1, std::vector<double> theta1,
                                                       int type1) {
        vec3<double> x = vec3<double>(pos1);
        vec3<double> th = vec3<double>(theta1);
        std::array<vec3<double>, 2> forceTorquex = forceTorque(x, th, type1);
        return msmrdtools::array2Dtovec2D(forceTorquex);
    }


    /* Evaluates force and torque. Template AUXVARIABLES corresponds to quaternion<double>, so it corresponds to
     * rigidbody particles with orientation described by a quaternion. */
    template<>
    std::vector<std::vector<double>>
    externalPotential<quaternion<double>>::forceTorquePyBind(std::vector<double> pos1, std::vector<double> theta1) {
        vec3<double> x = vec3<double>(pos1);
        quaternion<double> th = quaternion<double>(theta1);
        std::array<vec3<double>, 2> forceTorquex = forceTorque(x, th);
        return msmrdtools::array2Dtovec2D(forceTorquex);
    }


    /* Same as previous function with particle type dependence */
    template<>
    std::vector<std::vector<double>>
    externalPotential<quaternion<double>, int>::forceTorquePyBind(std::vector<double> pos1,
                                                                  std::vector<double> theta1, int type1) {
        vec3<double> x = vec3<double>(pos1);
        quaternion<double> th = quaternion<double>(theta1);
        std::array<vec3<double>, 2> forceTorquex = forceTorque(x, th, type1);
        return msmrdtools::array2Dtovec2D(forceTorquex);
    }


    /**
     *  Implementation of non-abstract functions of pairPotential abstract class (default constructor in header)
     *  All these functions are needed for PyBinding evaluate and forceTorque functions, since pyBind only accepts
     *  vectors as input. The different versions of the same functions depend on the AUXVARIABLES template that can
     *  be used for particles with no orientation, rod-like particles with their orientation described by one vectors,
     *  and by rigidbody particles with no symmetry and their orientation described by quaternions. AUXVARIABLES
     *  can also incorporate particle types to incorporate type-specific behavior in the potential function.
     */


    /* Evaluate functions */


    // Evaluates potential. Template AUXVARIABLES takes zero arguments, so it corresponds to no orientation.
    template<>
    double pairPotential<>::evaluatePyBind(std::vector<double> pos1, std::vector<double> pos2) {
        vec3<double> x1 = vec3<double>(pos1);
        vec3<double> x2 = vec3<double>(pos2);
        return evaluate(x1, x2);
    }

    /* Same as previous function with particle type dependence */
    template<>
    double pairPotential<int, int>::evaluatePyBind(std::vector<double> pos1, std::vector<double> pos2,
                                                   int type1, int type2) {
        vec3<double> x1 = vec3<double>(pos1);
        vec3<double> x2 = vec3<double>(pos2);
        return evaluate(x1, x2, type1, type2);
    }


    /* Evaluate potential. Template AUXVARIABLES corresponds to <vec3<double>, vec3<double>>, so it
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

    /* Same as previous function with particle type dependence */
    template<>
    double pairPotential<vec3<double>, vec3<double>, int, int>::evaluatePyBind(std::vector<double> pos1,
                                                                               std::vector<double> pos2,
                                                                               std::vector<double> theta1,
                                                                               std::vector<double> theta2,
                                                                               int type1, int type2) {
        vec3<double> x1 = vec3<double>(pos1);
        vec3<double> x2 = vec3<double>(pos2);
        vec3<double> th1 = vec3<double>(theta1);
        vec3<double> th2 = vec3<double>(theta2);
        return evaluate(x1, x2, th1, th2, type1, type2);
    }


    /* Evaluate potential. Template AUXVARIABLES corresponds to <quaternion<double>, quaternion<double>>, so it
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

    /* Same as previous function with particle type dependence */
    template<>
    double pairPotential<quaternion<double>, quaternion<double>, int, int>::evaluatePyBind(std::vector<double> pos1,
                                                                                           std::vector<double> pos2,
                                                                                           std::vector<double> theta1,
                                                                                           std::vector<double> theta2,
                                                                                           int type1, int type2) {
        vec3<double> x1 = vec3<double>(pos1);
        vec3<double> x2 = vec3<double>(pos2);
        quaternion<double> th1 = quaternion<double>(theta1);
        quaternion<double> th2 = quaternion<double>(theta2);
        return evaluate(x1, x2, th1, th2, type1, type2);
    }


    /* Force torque functions */


    // Evaluates force and torque. Template AUXVARIABLES takes zero arguments, so it corresponds to no orientation
    template<>
    std::vector<std::vector<double>> pairPotential<>::forceTorquePyBind(std::vector<double> pos1,
                                                                        std::vector<double> pos2) {
        vec3<double> x1 = vec3<double>(pos1);
        vec3<double> x2 = vec3<double>(pos2);
        std::array<vec3<double>, 4> forceTorquex = forceTorque(x1, x2);
        return msmrdtools::array2Dtovec2D(forceTorquex);
    }

    /* Same as previous function with particle type dependence */
    template<>
    std::vector<std::vector<double>> pairPotential<int, int>::forceTorquePyBind(std::vector<double> pos1,
                                                                                std::vector<double> pos2,
                                                                                int type1, int type2) {
        vec3<double> x1 = vec3<double>(pos1);
        vec3<double> x2 = vec3<double>(pos2);
        std::array<vec3<double>, 4> forceTorquex = forceTorque(x1, x2, type1, type2);
        return msmrdtools::array2Dtovec2D(forceTorquex);
    }


    /* Evaluates force and torque. Template AUXVARIABLES corresponds to <vec3<double>,vec3<double>>, so
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
        std::array<vec3<double>, 4> forceTorquex = forceTorque(x1, x2, th1, th2);
        return msmrdtools::array2Dtovec2D(forceTorquex);
    }

    /* Same as previous function with particle type dependence */
    template<>
    std::vector<std::vector<double>>
            pairPotential<vec3<double>, vec3<double>, int, int>::forceTorquePyBind(std::vector<double> pos1,
                                                                                   std::vector<double> pos2,
                                                                                   std::vector<double> theta1,
                                                                                   std::vector<double> theta2,
                                                                                   int type1, int type2) {
        vec3<double> x1 = vec3<double>(pos1);
        vec3<double> x2 = vec3<double>(pos2);
        vec3<double> th1 = vec3<double>(theta1);
        vec3<double> th2 = vec3<double>(theta2);
        std::array<vec3<double>, 4> forceTorquex = forceTorque(x1, x2, th1, th2, type1, type2);
        return msmrdtools::array2Dtovec2D(forceTorquex);
    }


    /* Evaluates force and torque. Template AUXVARIABLES corresponds to <quaternion<double>, quaternion<double>>,
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
        std::array<vec3<double>, 4> forceTorquex = forceTorque(x1, x2, th1, th2);
        return msmrdtools::array2Dtovec2D(forceTorquex);
    }

    /* Same as previous function with particle type dependence */
    template<>
    std::vector<std::vector<double>>
    pairPotential<quaternion<double>, quaternion<double>, int, int>::forceTorquePyBind(std::vector<double> pos1,
                                                                                       std::vector<double> pos2,
                                                                                       std::vector<double> theta1,
                                                                                       std::vector<double> theta2,
                                                                                       int type1, int type2) {
        vec3<double> x1 = vec3<double>(pos1);
        vec3<double> x2 = vec3<double>(pos2);
        quaternion<double> th1 = quaternion<double>(theta1);
        quaternion<double> th2 = quaternion<double>(theta2);
        std::array<vec3<double>, 4> forceTorquex = forceTorque(x1, x2, th1, th2, type1, type2);
        return msmrdtools::array2Dtovec2D(forceTorquex);
    }

}