//
// Created by maojrs on 9/26/18.
//

#include "potentials/patchyProtein.hpp"
#include "quaternion.hpp"


namespace msmrd {
   /*
    * Constructors analogous to the patchyParticle potential.
    * @param sigma diameter of sphere at which patches are placed. Corresponds to the radius of the particle.
    * @param patchesCoordintates list of patches coordinates in a sphere of unit radius
    */
    patchyProtein::patchyProtein(double sigma1, double sigma2, double strength1, double strength2,
                                 std::vector<vec3<double>> patchesCoordinates)
            :  sigma1(sigma1), sigma2(sigma2), strength1(strength1), strength2(strength2),
               patchesCoordinates(patchesCoordinates) {
       setPotentialParameters();
    };

    patchyProtein::patchyProtein(double sigma1, double sigma2, double strength1, double strength2,
                                 std::vector<std::vector<double>> patchesCoords)
            : sigma1(sigma1), sigma2(sigma2), strength1(strength1), strength2(strength2) {
        patchesCoordinates.resize(patchesCoords.size());
        for (int i=0; i < patchesCoords.size(); i++ ) {
            patchesCoordinates[i] = vec3<double> (patchesCoords[i]);
        }
        setPotentialParameters();
    }


    // Set potential parameters for potential
    void patchyProtein::setPotentialParameters() {
        // Set strengths of the three potential parts
        epsRepulsive[0] = 1.0*strength1;
        epsAttractive[0] = -0.0*strength1;
        epsPatches[0] = -0.15*strength1;

        // Set stiffness for the three potentials
        aRepulsive[0] = 1.5;
        aAttractive[0] = 0.75;
        aPatches[0] = 40.0;

        // Set range parameter for the three potentials
        rstarRepulsive[0] = 0.75*sigma1;
        rstarAttractive[0] = 0.85*sigma1;
        rstarPatches[0] = 0.1*sigma1;

        // Set strengths of the three potential parts
        epsRepulsive[1] = 1.0*strength2;
        epsAttractive[1] = -0.0*strength2;
        epsPatches[1] = -0.15*strength2;

        // Set stiffness for the three potentials
        aRepulsive[1] = 1.5;
        aAttractive[1] = 0.75;
        aPatches[1] = 40.0;

        // Set range parameter for the three potentials
        rstarRepulsive[1] = 0.75*sigma2;
        rstarAttractive[1] = 0.85*sigma2;
        rstarPatches[1] = 0.1*sigma2;
    }


    // Evaluates potential at given positions and orientations of two particles
    double patchyProtein::evaluate(vec3<double> pos1, vec3<double> pos2,
                                   quaternion<double> theta1, quaternion<double> theta2,
                                    int type1, int type2) {
        double repulsivePotential;
        double attractivePotential;
        double patchesPotential = 0.0;
        vec3<double> patch1;
        vec3<double> patch2;
        vec3<double> patchNormal1;
        vec3<double> patchNormal2;
        vec3<double> rpatch;
        vec3<double> rvec = pos2 - pos1;

        repulsivePotential = quadraticPotential(rvec.norm(), sigma1, epsRepulsive[0], aRepulsive[0], rstarRepulsive[0]);
        attractivePotential = quadraticPotential(rvec.norm(), sigma1, epsAttractive[0], aAttractive[0], rstarAttractive[0]);

        if (rvec.norm() <= 2*sigma1) {
            // Loop over all patches
            for (int i = 0; i < patchesCoordinates.size(); i++) {
                patchNormal1 = rotateVec(patchesCoordinates[i], theta1);
                patch1 = pos1 + 0.5*sigma1*patchNormal1;
                for (int j = 0; j < patchesCoordinates.size(); j++) {
                    patchNormal2 = rotateVec(patchesCoordinates[j], theta2);
                    patch2 = pos2 + 0.5*sigma1*patchNormal2;
                    rpatch = patch2 - patch1; // Scale unit distance of patches by sigma
                    patchesPotential += quadraticPotential(rpatch.norm(), sigma1, epsPatches[0], aPatches[0], rstarPatches[0]);
                }
            }
        }
        return repulsivePotential + attractivePotential + patchesPotential;
    }

    /* Calculate and return (force1, torque1, force2, torque2), which correspond to the force and torque
     * acting on particle1 and the force and torque acting on particle2, respectively. */
    std::array<vec3<double>, 4> patchyProtein::forceTorque(vec3<double> pos1, vec3<double> pos2,
                                                           quaternion<double> theta1, quaternion<double> theta2,
                                                           int type1, int type2) {
        vec3<double> force;
        vec3<double> force1 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> force2 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> torque1 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> torque2 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> rvec = pos2 - pos1;
        // auxiliary variables to calculate force and torque
        double repulsiveForceNorm;
        double attractiveForceNorm;
        double patchesForceNorm;
        vec3<double> patchForce;
        vec3<double> patch1;
        vec3<double> patch2;
        vec3<double> patchNormal1;
        vec3<double> patchNormal2;
        vec3<double> rpatch;

        /* Calculate and add forces due to repulsive and attractive isotropic potentials.
         *  Note correct sign/direction of force given by rvec/rvec.norm*() */
        repulsiveForceNorm = derivativeQuadraticPotential(rvec.norm(), sigma1, epsRepulsive[0], aRepulsive[0], rstarRepulsive[0]);
        attractiveForceNorm = derivativeQuadraticPotential(rvec.norm(), sigma1, epsAttractive[0], aAttractive[0], rstarAttractive[0]);
        force = (repulsiveForceNorm + attractiveForceNorm)*rvec/rvec.norm();

        // Calculate forces and torque due to patches interaction
        if (rvec.norm() <= 2*sigma1) {
            // Loop over all patches of particle 1
            for (int i = 0; i < patchesCoordinates.size(); i++) {
                patchNormal1 = rotateVec(patchesCoordinates[i], theta1);
                patchNormal1 = patchNormal1/patchNormal1.norm();
                patch1 = pos1 + 0.5*sigma1*patchNormal1;
                // Loop over all patches of particle 2
                for (int j = 0; j < patchesCoordinates.size(); j++) {
                    patchNormal2 = rotateVec(patchesCoordinates[j], theta2);
                    patchNormal2 = patchNormal2/patchNormal2.norm();
                    patch2 = pos2 + 0.5*sigma1*patchNormal2;
                    rpatch = patch2 - patch1;
                    // Calculate force vector between patches , correct sign of force given by rpatch/rpatch.norm().
                    patchesForceNorm = derivativeQuadraticPotential(rpatch.norm(), sigma1, epsPatches[0], aPatches[0], rstarPatches[0]);
                    if (rpatch.norm() == 0) {
                        patchForce = vec3<double> (0, 0, 0);
                    }
                    else {
                        patchForce = patchesForceNorm*rpatch/rpatch.norm();
                    }
                    // Calculate force and torque acting on particle 1 and add values to previous forces and torques
                    force1 += patchForce;
                    torque1 += 0.5*sigma1 * patchNormal1.cross(patchForce);

                    // Calculate force and torque acting on particle 2 and add values to previous forces and torques
                    force2 += -1.0*patchForce;
                    torque2 += 0.5*sigma1 * patchNormal2.cross(-1.0*patchForce);
                }
            }
        }
        return {force + force1, torque1, -1.0*force + force2, torque2};
    }


    /* Custom quadratic potential (same functions as in patchy particle, inheritance was not an option :( ):
     * @param r is distance between particles or  patches,
     * @param eps strength of the potential
     * @param rstar distance to switch to second potential function
     * @param a the stiffness, together with rstar determines the range of the potential */
    double patchyProtein::quadraticPotential(double r, double sig, double eps, double a, double rstar) {
        // Parameters to force continuity and continuous derivative
        double rcritical = std::pow(sig, 2)/(a*rstar);
        double b = (1.0 - a*std::pow(rstar/sig, 2))/std::pow(rcritical/sig - rstar/sig, 2);
        // Distance cannot be negative, so it is always the absolute value
        double rlocal = std::abs(r);
        // Calculate different regions of potential
        if (rlocal < rstar) {
            return eps * (1.0 - a * std::pow((rlocal / sig),2));
        } else if (rlocal < rcritical) {
            return eps * b * std::pow((rcritical / sig - rlocal / sig), 2);
        } else {
            return 0.0;
        }
    }

    // Derivative of patchyProtein::quadraticPotential
    double patchyProtein::derivativeQuadraticPotential(double r, double sig, double eps, double a, double rstar) {
        // Parameters to force continuity and continuous derivative
        double rcritical = std::pow(sig, 2)/(a*rstar);
        double b = (1.0 - a*std::pow(rstar/sig, 2))/std::pow(rcritical/sig - rstar/sig, 2);
        // Distance cannot be negative, so it is always the absolute value
        double rlocal = std::abs(r);
        // Calculate derivative in different regions of potential
        if (rlocal < rstar) {
            return - 2.0*eps*a * rlocal/std::pow(sig,2);
        } else if (rlocal < rcritical) {
            return -2.0*eps * b * (rcritical - rlocal)/std::pow(sig,2);
        } else {
            return 0.0;
        }
    }


}