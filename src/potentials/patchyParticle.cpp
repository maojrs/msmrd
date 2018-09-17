//
// Created by maojrs on 9/10/18.
//

#include "potentials/patchyParticle.hpp"
#include "quaternion.hpp"




namespace msmrd{ ;
    /*
     * @param radius radius at which patches are placed. Corresponds to the radius of the particle
     * @param patchesCoordintates list of patches coordinates in a sphere of unit radius
     */
    patchyParticle::patchyParticle(double sigma, double strength, std::vector<vec3<double>> patchesCoordinates)
            : sigma(sigma), strength(strength), patchesCoordinates(patchesCoordinates)  {
        // Set strengths of the three potential parts
        epsRepulsive = strength;
        epsPatches = -0.2*strength;
        epsAttractive = -0.1*strength;

        // Set stiffness for the three potentials
        aRepulsive = 1.0;
        aPatches = 20.0;
        aAttractive = 0.5;

        // Set range parameter for the three potentials
        rstarRepulsive = 0.85*sigma;
        rstarPatches = 0.1*sigma;
        rstarAttractive = 0.85*sigma;
    }

    patchyParticle::patchyParticle(double sigma, double strength, std::vector<std::vector<double>> patchesCoords)
            : sigma(sigma), strength(strength) {
        patchesCoordinates.resize(patchesCoords.size());
        for (int i=0; i < patchesCoords.size(); i++ ) {
            patchesCoordinates[i] = vec3<double> (patchesCoords[i]);
        }
        // Set strengths of the three potential parts
        epsRepulsive = strength;
        epsAttractive = -0.1*strength;
        epsPatches = -0.2*strength;

        // Set stiffness for the three potentials
        aRepulsive = 1.0;
        aAttractive = 0.5;
        aPatches = 20.0;

        // Set range parameter for the three potentials
        rstarRepulsive = 0.85*sigma;
        rstarAttractive = 0.85*sigma;
        rstarPatches = 0.1*sigma;
    }


    double patchyParticle::evaluate(vec3<double> pos1, vec3<double> pos2, quaternion<double> theta1, quaternion<double> theta2) {
        double repulsivePotential = 0.0;
        double attractivePotential = 0.0;
        double patchesPotential = 0.0;
        vec3<double> patch1;
        vec3<double> patch2;
        vec3<double> rpatch;
        vec3<double> rvec = pos2 - pos1;

        repulsivePotential = quadraticPotential(rvec.norm(), epsRepulsive, aRepulsive, rstarRepulsive);
        attractivePotential = quadraticPotential(rvec.norm(), epsAttractive, aAttractive, rstarAttractive);

        if (rvec.norm() <= 2*sigma) {
            // Loop over all patches
            for (int i = 0; i < patchesCoordinates.size(); i++) {
                patch1 = rotateVec(patchesCoordinates[i], theta1);
                patch1 = pos1 + 0.5*sigma*patch1;
                for (int j = 0; j < patchesCoordinates.size(); j++) {
                    patch2 = rotateVec(patchesCoordinates[j], theta2);
                    patch2 = pos2 + 0.5*sigma*patch2;
                    rpatch = patch1 - patch2; // Scale unit distance of patches by sigma
                    patchesPotential += quadraticPotential(rpatch.norm(), epsPatches, aPatches, rstarPatches);
                }
            }
        }
        return repulsivePotential + attractivePotential + patchesPotential;
        //return patchesPotential;
    }

    std::array<vec3<double>, 2> patchyParticle::forceTorque(vec3<double> pos1, vec3<double> pos2, quaternion<double> theta1, quaternion<double> theta2) {
        vec3<double> force = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> torque = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> rvec = pos2 - pos1;
        // auxiliary variables to calculate force and torque
        double repulsiveForceNorm = 0.0;
        double attractiveForceNorm = 0.0;
        double patchesForceNorm = 0.0;
        vec3<double> patchForce;
        vec3<double> patchForceNormal;
        vec3<double> patchForceTangent;
        vec3<double> patch1;
        vec3<double> patch2;
        vec3<double> patchNormal1;
        vec3<double> patchNormal2;
        vec3<double> rpatch;

        // Calculate and add forces due to repulsive and attractive isotropic potentials.
        repulsiveForceNorm = -1.0*derivativeQuadraticPotential(rvec.norm(), epsRepulsive, aRepulsive, rstarRepulsive);
        attractiveForceNorm = -1.0*derivativeQuadraticPotential(rvec.norm(), epsAttractive, aAttractive, rstarAttractive);
        force = (repulsiveForceNorm + attractiveForceNorm)*rvec/rvec.norm();

        // Calculate forces and torque due to patches interaction
        if (rvec.norm() <= 2*sigma) {
            // Loop over all patches of particle 1
            for (int i = 0; i < patchesCoordinates.size(); i++) {
                patchNormal1 = rotateVec(patchesCoordinates[i], theta1);
                patch1 = pos1 + 0.5*sigma*patchNormal1;
                // Loop over all patches of particle 2
                for (int j = 0; j < patchesCoordinates.size(); j++) {
                    patchNormal2 = rotateVec(patchesCoordinates[j], theta2);
                    patch2 = pos2 + 0.5*sigma*patchNormal2;
                    rpatch = patch2 - patch1; // Scale unit distance of patches by sigma
                    // Calculate force vector between patches
                    patchesForceNorm = -1.0*derivativeQuadraticPotential(rpatch.norm(), epsPatches, aPatches, rstarPatches);
                    patchForce = patchesForceNorm*rpatch/rpatch.norm();
                    // Decompose into tangential and normal components according to normal at patch
                    patchForceNormal = (patchForce * patchNormal1) * patchNormal1;
                    patchForceTangent = patchForce - patchForceNormal;
                    // Add values to previous forces and torques
                    force += patchForceNormal;
                    torque += 0.5*sigma * patchNormal1.cross(patchForceTangent);
                }
            }
        }
        return {force, torque};
    }

    /* Custom quadratic potential:
     * @param r is distance between particles or  patches,
     * @param eps strength of the potential
     * @param rstar distance to switch to second potential function
     * @param a the stiffness, together with rstar determines the range of the potential */
    double patchyParticle::quadraticPotential(double r, double eps, double a, double rstar) {
        // Parameters to force continuity and continuous derivative
        double rcritical = std::pow(sigma, 2)/(a*rstar);
        double b = (1.0 - a*std::pow(rstar/sigma, 2))/std::pow(rcritical/sigma - rstar/sigma, 2);
        // Distance cannot be negative, so it is always the absolute value
        double rlocal = std::abs(r);
        // Calculate different regions of potential
        if (rlocal < rstar) {
            return eps * (1.0 - a * std::pow((rlocal / sigma),2));
        } else if (rlocal < rcritical) {
            return eps * b * std::pow((rcritical / sigma - rlocal / sigma), 2);
        } else {
            return 0.0;
        }
    }

    double patchyParticle::derivativeQuadraticPotential(double r, double eps, double a, double rstar) {
        // Parameters to force continuity and continuous derivative
        double rcritical = std::pow(sigma, 2)/(a*rstar);
        double b = (1.0 - a*std::pow(rstar/sigma, 2))/std::pow(rcritical/sigma - rstar/sigma, 2);
        // Distance cannot be negative, so it is always the absolute value
        double rlocal = std::abs(r);
        // Calculate derivative in different regions of potential
        if (r < rstar) {
            return - 2.0*eps*a * r/std::pow(sigma,2);
        } else if (r < rcritical) {
            return -2.0*eps * b * (rcritical - r)/std::pow(sigma,2);
        } else {
            return 0.0;
        }
    }

}

