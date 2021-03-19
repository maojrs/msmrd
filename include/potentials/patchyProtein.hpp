//
// Created by maojrs on 9/26/18.
//

#pragma once
#include "potentials.hpp"

namespace msmrd{

    /*
     * Declaration of potential function for patchy protein model. This class is a child from the
     * general pair potentials class that accepts quaternion-based orientation (rigidBody particle type).
     * Some of its functionality is based on the patchy particles class. The main difference is that this class
     * can take different type of patches. Note the potential will depend on the position of both particles,
     * their orientation (quaternion<double>, quaternion<double), and their types (int, int).
     */
    class patchyProtein : public pairPotential {
    private:
        // Need to call this function from constructor to define parameters
        void setPotentialParameters();

    protected:
        std::vector<vec3<double>> patchesCoordinatesA;
        std::vector<vec3<double>> patchesCoordinatesB;
        bool patchesActive = true;

        double sigma = 1.0;
        double strength = 100.0;
        double patchPotentialScaling = 1.0;

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

        double quadraticPotential(double r, double sig, double eps, double a, double rstar);

        double derivativeQuadraticPotential(double r, double sig, double eps, double a, double rstar);

        virtual std::vector<vec3<double>> assignPatches(int type);

        virtual double evaluatePatchesPotential( particle &part1, particle &part2,
                                        vec3<double> &pos1virtual,
                                        std::vector<vec3<double>> &patchesCoords1,
                                        std::vector<vec3<double>> &patchesCoords2);

        virtual std::array<vec3<double>, 4> forceTorquePatches(particle &part1, particle &part2,
                                                       vec3<double> &pos1virtual,
                                                       std::vector<vec3<double>> &patchesCoords1,
                                                       std::vector<vec3<double>> &patchesCoords2);

    public:
        patchyProtein() = default;

        patchyProtein(double sigma, double strength,
                      std::vector<vec3<double>> patchesCoordinatesA,
                      std::vector<vec3<double>> patchesCoordinatesB);
        patchyProtein(double sigma, double strength,
                      std::vector<std::vector<double>> patchesCoordinatesA,
                      std::vector<std::vector<double>> patchesCoordinatesB);

        double evaluate(particle &part1, particle &part2) override;

        std::array<vec3<double>, 4> forceTorque(particle &part1, particle &part2) override;

        bool arePatchesActive() { return patchesActive; }

    };


}

