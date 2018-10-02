//
// Created by maojrs on 9/10/18.
//

#include "potentials/patchyParticle.hpp"
#include "quaternion.hpp"
#include "tools.hpp"




namespace msmrd{ ;
    /*
     * @param sigma diameter of sphere at which patches are placed. Corresponds to the radius of the particle.
     * @param strength overall strength of the potential
     * @param patchesCoordinates list of patches coordinates in a sphere of unit radius
     */
    patchyParticle::patchyParticle(double sigma, double strength, std::vector<vec3<double>> patchesCoordinates)
            : sigma(sigma), strength(strength), patchesCoordinates(patchesCoordinates)  {

        // Set strengths of the three potential parts
        epsRepulsive = 1.0*strength;
        epsAttractive = -0.0*strength;
        epsPatches = -0.15*strength;

        // Set stiffness for the three potentials
        aRepulsive = 1.5;
        aAttractive = 0.75;
        aPatches = 40.0;

        // Set range parameter for the three potentials
        rstarRepulsive = 0.75*sigma;
        rstarAttractive = 0.85*sigma;
        rstarPatches = 0.1*sigma;

        for (auto &patch : patchesCoordinates) {
            if (patch.norm() != 1) {
                throw std::range_error("Patches coordinates must be in the unit sphere");
            }
        }
    }

    patchyParticle::patchyParticle(double sigma, double strength, std::vector<std::vector<double>> patchesCoords)
            : sigma(sigma), strength(strength) {
        patchesCoordinates.resize(patchesCoords.size());
        for (int i=0; i < patchesCoords.size(); i++ ) {
            patchesCoordinates[i] = vec3<double> (patchesCoords[i]);
        }

        // Set strengths of the three potential parts
        epsRepulsive = 1.0*strength;
        epsAttractive = -0.0*strength;
        epsPatches = -0.3*strength;

        // Set stiffness for the three potentials
        aRepulsive = 1.5;
        aAttractive = 0.75;
        aPatches = 40.0;

        // Set range parameter for the three potentials
        rstarRepulsive = 0.75*sigma;
        rstarAttractive = 0.85*sigma;
        rstarPatches = 0.1*sigma;

        for (auto &patch : patchesCoordinates) {
            if (patch.norm() != 1) {
                throw std::range_error("Patches coordinates must be in the unit sphere");
            }
        }
    }

    // Evaluates potential at given positions and orientations of two particles
    double patchyParticle::evaluate(vec3<double> pos1, vec3<double> pos2, quaternion<double> theta1, quaternion<double> theta2) {
        double repulsivePotential;
        double attractivePotential;
        double patchesPotential = 0.0;
        double angleModulation;
        vec3<double> patch1;
        vec3<double> patch2;
        vec3<double> patchNormal1;
        vec3<double> patchNormal2;
        vec3<double> rpatch;
        vec3<double> rvec = pos2 - pos1;

        repulsivePotential = quadraticPotential(rvec.norm(), sigma, epsRepulsive, aRepulsive, rstarRepulsive);
        attractivePotential = quadraticPotential(rvec.norm(), sigma, epsAttractive, aAttractive, rstarAttractive);

        if (rvec.norm() <= 2*sigma) {
            // Loop over all patches
            for (int i = 0; i < patchesCoordinates.size(); i++) {
                patchNormal1 = msmrdtools::rotateVec(patchesCoordinates[i], theta1);
                patch1 = pos1 + 0.5*sigma*patchNormal1;
                for (int j = 0; j < patchesCoordinates.size(); j++) {
                    patchNormal2 = msmrdtools::rotateVec(patchesCoordinates[j], theta2);
                    patch2 = pos2 + 0.5*sigma*patchNormal2;
                    rpatch = patch2 - patch1; // Scale unit distance of patches by sigma
//                    angleModulation = -1*patchNormal1*patchNormal2; //angular modulation
//                    if (angleModulation < 0.0) {
//                        angleModulation = 0.0;
//                    }
//                    angleModulation*
                    patchesPotential += quadraticPotential(rpatch.norm(), sigma, epsPatches, aPatches, rstarPatches);
                }
            }
        }
        return repulsivePotential + attractivePotential + patchesPotential;
    }

    /* Calculate and return (force1, torque1, force2, torque2), which correspond to the force and torque
     * acting on particle1 and the force and torque acting on particle2, respectively. */
    std::array<vec3<double>, 4> patchyParticle::forceTorque(vec3<double> pos1, vec3<double> pos2, quaternion<double> theta1, quaternion<double> theta2) {
        vec3<double> force;
        vec3<double> force1 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> force2 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> torque1 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> torque2 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> rvec = pos2 - pos1;
        // auxiliary variables to calculate force and torque
        double repulsiveForceNorm = 0.0;
        double attractiveForceNorm = 0.0;
        double patchesForceNorm = 0.0;
        double angleModulation;
        vec3<double> patchForce;
        vec3<double> patchTorque;
        vec3<double> patch1;
        vec3<double> patch2;
        vec3<double> patchNormal1;
        vec3<double> patchNormal2;
        vec3<double> rpatch;

        /* Calculate and add forces due to repulsive and attractive isotropic potentials.
         *  Note correct sign/direction of force given by rvec/rvec.norm*() */
        repulsiveForceNorm = derivativeQuadraticPotential(rvec.norm(), sigma, epsRepulsive, aRepulsive, rstarRepulsive);
        attractiveForceNorm = derivativeQuadraticPotential(rvec.norm(), sigma, epsAttractive, aAttractive, rstarAttractive);
        force = (repulsiveForceNorm + attractiveForceNorm)*rvec/rvec.norm();

        // Calculate forces and torque due to patches interaction
        if (rvec.norm() <= 2*sigma) {
            // Loop over all patches of particle 1
            for (int i = 0; i < patchesCoordinates.size(); i++) {
                patchNormal1 = msmrdtools::rotateVec(patchesCoordinates[i], theta1);
                patchNormal1 = patchNormal1/patchNormal1.norm();
                patch1 = pos1 + 0.5*sigma*patchNormal1;
                // Loop over all patches of particle 2
                for (int j = 0; j < patchesCoordinates.size(); j++) {
                    patchNormal2 = msmrdtools::rotateVec(patchesCoordinates[j], theta2);
                    patchNormal2 = patchNormal2/patchNormal2.norm();
                    patch2 = pos2 + 0.5*sigma*patchNormal2;
                    rpatch = patch2 - patch1;
                    // Calculate force vector between patches , correct sign of force given by rpatch/rpatch.norm().
                    patchesForceNorm = derivativeQuadraticPotential(rpatch.norm(), sigma, epsPatches, aPatches, rstarPatches);
                    if (rpatch.norm() == 0) {
                        patchForce = vec3<double> (0, 0, 0);
                    }
                    else {
//                        angleModulation = -1*patchNormal1*patchNormal2;
//                        if (angleModulation < 0.0) {
//                            angleModulation = 0.0;
//                        }
//                        patchesForceNorm = angleModulation*patchesForceNorm;
                        patchForce = patchesForceNorm*rpatch/rpatch.norm();
                    }
                    // Calculate force and torque acting on particle 1 and add values to previous forces and torques
                    force1 += patchForce;
                    torque1 += 0.5*sigma * patchNormal1.cross(patchForce);

//                    // Caluclate additional term from angleModulator product rule
//                    patchTorque = patchNormal1.cross(patchNormal2);
//                    patchTorque = patchTorque*quadraticPotential(rpatch.norm(), sigma, epsPatches, aPatches, rstarPatches);
//                    torque += patchTorque;

                    // Calculate force and torque acting on particle 2 and add values to previous forces and torques
                    force2 += -1.0*patchForce;
                    torque2 += 0.5*sigma * patchNormal2.cross(-1.0*patchForce);
                }
            }
        }
        return {force + force1, torque1, -1.0*force + force2, torque2};
    }

    /* Custom quadratic potential:
     * @param r is distance between particles or  patches,
     * @param eps strength of the potential
     * @param rstar distance to switch to second potential function
     * @param a the stiffness, together with rstar determines the range of the potential */
    double patchyParticle::quadraticPotential(double r, double sig, double eps, double a, double rstar) {
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

    // Derivative of patchyParticle::quadraticPotential
    double patchyParticle::derivativeQuadraticPotential(double r, double sig, double eps, double a, double rstar) {
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

