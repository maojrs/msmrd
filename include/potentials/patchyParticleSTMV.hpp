//
// Created by maojrs on 4/29/21.
//
#pragma once
#include "potentials/patchyParticleAngular.hpp"
#include "quaternion.hpp"
#include "tools.hpp"

namespace msmrd {
    /*
     * Declaration of potential function for modeling the Satellite Tobacco Mosaic Virus (STMV).
     * It is based on patchy particles approach. The repulsive potential function avoids
     * overlapping of a complex shape made of three spherical parts (m1,m2,m3). The attraction
     * sites or patches are located around the surface of this three spheres and are denoted by
     * i1 to i5. A brief description of the structur follows:
     * In its ground orientation (values in angstroms (10^(-10)m)):
     * m1 position: [18.15254758, 0, 0] radius: 8.56
     * m2 position: [0, 0, 0] radius: 13.58  **reference position**
     * m3 position: [-5.86991122, 14.60740195, 0] radius: 8.18
     * i1 position: [16.94500527, 6.74187053, -6.94570924]
     * i2 position: [21.98599972, -2.07424062, -8.26036629]
     * i3 position: [-4.87943602, -9.24886815, -1.63915576]
     * i4 position: [-13.48049652, 8.41422688, -0.57270532]
     * i5 position: [-10.03689633, 17.39415955, -5.62472092]
     * The radii of all patches is 4, but they produce no non-overlapping potential.
     * Note the implementation will be in nanometers (10^(-9)m)
     * There are three possible type of interaction of two identical particles:
     * interaction 1: i1:i2
     * interaction 2: i3:i4
     * interaction 5: i5:i5
     */
    class patchyParticleSTMV : public pairPotential {
    private:
        // Need to call this function from constructor to define parameters
        void setPotentialParameters();

    protected:
        std::vector<vec3<double>> structureCoordinates;
        std::vector<vec3<double>> patchesCoordinates;

        std::array<double, 3> sigma;
        double sigmaPatches;
        double strength = 100.0;
        double angularStrength = 2.0;
        double patchPotentialScaling = 1.0;

        // Parameters for first and second type of interactions
        // Strength of potentials
        double epsRepulsive;
        std::array<double, 3> epsPatches;
        // Stiffness of potentials
        double aRepulsive;
        std::array<double, 3> aPatches;
        // Range parameter of potentials
        double rstarRepulsive;
        std::array<double, 3> rstarPatches;

        double quadraticPotential(double r, double sig, double eps, double a, double rstar);

        double derivativeQuadraticPotential(double r, double sig, double eps, double a, double rstar);

        std::array<vec3<double>, 4> forceTorquePatches(
                particle &part1, particle &part2, const vec3<double> pos1virtual);

        std::array<vec3<double>, 4> forceTorqueStructure(
                particle &part1, particle &part2, const vec3<double> pos1virtual);

    public:
        patchyParticleSTMV(double strength, double angularStrength);

        double evaluate(particle &part1, particle &part2) override;

        std::array<vec3<double>, 4>
        forceTorque(particle &part1, particle &part2) override;

        std::array<double, 3> getPartPosition(std::string particlePart, particle &part);

    };



}
