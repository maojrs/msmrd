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
        std::vector<vec3<double>> patchesCoordinatesA;
        std::vector<vec3<double>> patchesCoordinatesB;
        double sigma = 1.0;
        double strength = 100.0;

        // Parameters for first and second type of interactions
        // Strength of potentials
        double epsRepulsive;
        double epsAttractive;
        std::array<double, 2> epsPatches;
        // Stiffness of potentials
        double aRepulsive;
        double aAttractive;
        std::array<double, 2> aPatches;
        // Range parameter of potentials
        double rstarRepulsive;
        double rstarAttractive;
        std::array<double, 2> rstarPatches;

        void setPotentialParameters();
        std::vector<vec3<double>> assignPatches(int type);
        double quadraticPotential(double r, double sig, double eps, double a, double rstar);
        double derivativeQuadraticPotential(double r, double sig, double eps, double a, double rstar);

    public:
        patchyProtein(double sigma, double strength,
                      std::vector<vec3<double>> patchesCoordinatesA,
                      std::vector<vec3<double>> patchesCoordinatesB);
        patchyProtein(double sigma, double strength,
                      std::vector<std::vector<double>> patchesCoordinatesA,
                      std::vector<std::vector<double>> patchesCoordinatesB);

        double evaluate(vec3<double> pos1, vec3<double> pos2,
                        quaternion<double> theta1, quaternion<double> theta2,
                        int type1, int type2) override;

        std::array<vec3<double>, 4> forceTorque(vec3<double> pos1, vec3<double> pos2,
                                                quaternion<double> theta1, quaternion<double> theta2,
                                                int type1, int type2) override;

    };


}
