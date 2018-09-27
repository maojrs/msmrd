//
// Created by maojrs on 9/26/18.
//

#pragma once
#include "potentials/patchyParticle.hpp"

namespace msmrd{

    /*
     * Declaration of potential function for patchy protein model. This class is a child from the
     * general pair potentials class that accepts quaternion-based orientation and particle types
     * (rigidbodyMix). Some of its functionality is based on the patchy particles class. Note the
     * potential will depend on the position of both particles, their orientation (quaternion<double>,
     * quaternion<double), and their types (int, int).
     */
    class patchyProtein : public pairPotential<quaternion<double>, quaternion<double>, int, int> {
    private:
        std::vector<vec3<double>> patchesCoordinates;
        double sigma1 = 1.0;
        double sigma2 = 1.0;
        double strength1 = 100.0;
        double strength2 = 100.0;

        // Parameters for first and second type of interactions
        // Strength of potentials
        std::array<double, 2> epsRepulsive;
        std::array<double, 2> epsAttractive;
        std::array<double, 2> epsPatches;
        // Stiffness of potentials
        std::array<double, 2> aRepulsive;
        std::array<double, 2> aAttractive;
        std::array<double, 2> aPatches;
        // Range parameter of potentials
        std::array<double, 2> rstarRepulsive;
        std::array<double, 2> rstarAttractive;
        std::array<double, 2> rstarPatches;

        void setPotentialParameters();
        double quadraticPotential(double r, double sig, double eps, double a, double rstar);
        double derivativeQuadraticPotential(double r, double sig, double eps, double a, double rstar);

    public:
        patchyProtein(double sigma1, double sigma2, double strength1, double strength2,
                      std::vector<vec3<double>> patchesCoordinates);
        patchyProtein(double sigma1, double sigma2, double strength1, double strength2,
                      std::vector<std::vector<double>> patchesCoordinates);

        double evaluate(vec3<double> pos1, vec3<double> pos2,
                        quaternion<double> theta1, quaternion<double> theta2,
                        int type1, int type2) override;

        std::array<vec3<double>, 4>
        forceTorque(vec3<double> pos1, vec3<double> pos2,
                    quaternion<double> theta1, quaternion<double> theta2,
                    int type1, int type2) override;

    };


}

