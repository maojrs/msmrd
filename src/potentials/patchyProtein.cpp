//
// Created by maojrs on 9/26/18.
//

#include "potentials/patchyProtein.hpp"
#include "quaternion.hpp"


namespace msmrd {
   /*
    * Constructors analogous to the patchyParticle potential.
    * @param radius radius at which patches are placed. Corresponds to the radius of the particle
    * @param patchesCoordintates list of patches coordinates in a sphere of unit radius
    */
    patchyProtein::patchyProtein(double sigma, double strength, std::vector<vec3<double>> patchesCoordinates)
            : patchyParticle(sigma, strength, patchesCoordinates) {
        // Set patchy particle parameters (strength, stiffness and range) for attractive isotropic potential
        epsAttractive = -0.05*strength;
        aAttractive = 0.75;
        rstarAttractive = 0.85*sigma;
    };

    patchyProtein::patchyProtein(double sigma, double strength, std::vector<std::vector<double>> patchesCoordinates)
            : patchyParticle(sigma, strength, patchesCoordinates) {
        // Set patchy particle parameters (strength, stiffness and range) for attractive isotropic potential
        epsAttractive = -0.05*strength;
        aAttractive = 0.75;
        rstarAttractive = 0.85*sigma;
    }

        // Evaluates potential at given positions and orientations of two particles
    double patchyProtein::evaluate(vec3<double> pos1, vec3<double> pos2, quaternion<double> theta1, quaternion<double> theta2) {
        double repulsivePotential;
        double attractivePotential;
        double patchesPotential = 0.0;
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
                patchNormal1 = rotateVec(patchesCoordinates[i], theta1);
                patch1 = pos1 + 0.5*sigma*patchNormal1;
                for (int j = 0; j < patchesCoordinates.size(); j++) {
                    patchNormal2 = rotateVec(patchesCoordinates[j], theta2);
                    patch2 = pos2 + 0.5*sigma*patchNormal2;
                    rpatch = patch2 - patch1; // Scale unit distance of patches by sigma
                    patchesPotential += quadraticPotential(rpatch.norm(), sigma, epsPatches, aPatches, rstarPatches);
                }
            }
        }
        return repulsivePotential + attractivePotential + patchesPotential;
    }

    /* Calculate and return (force1, torque1, force2, torque2), which correspond to the force and torque
     * acting on particle1 and the force and torque acting on particle2, respectively. */
    std::array<vec3<double>, 4> patchyProtein::forceTorque(vec3<double> pos1, vec3<double> pos2, quaternion<double> theta1, quaternion<double> theta2) {
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
        repulsiveForceNorm = derivativeQuadraticPotential(rvec.norm(), sigma, epsRepulsive, aRepulsive, rstarRepulsive);
        attractiveForceNorm = derivativeQuadraticPotential(rvec.norm(), sigma, epsAttractive, aAttractive, rstarAttractive);
        force = (repulsiveForceNorm + attractiveForceNorm)*rvec/rvec.norm();

        // Calculate forces and torque due to patches interaction
        if (rvec.norm() <= 2*sigma) {
            // Loop over all patches of particle 1
            for (int i = 0; i < patchesCoordinates.size(); i++) {
                patchNormal1 = rotateVec(patchesCoordinates[i], theta1);
                patchNormal1 = patchNormal1/patchNormal1.norm();
                patch1 = pos1 + 0.5*sigma*patchNormal1;
                // Loop over all patches of particle 2
                for (int j = 0; j < patchesCoordinates.size(); j++) {
                    patchNormal2 = rotateVec(patchesCoordinates[j], theta2);
                    patchNormal2 = patchNormal2/patchNormal2.norm();
                    patch2 = pos2 + 0.5*sigma*patchNormal2;
                    rpatch = patch2 - patch1;
                    // Calculate force vector between patches , correct sign of force given by rpatch/rpatch.norm().
                    patchesForceNorm = derivativeQuadraticPotential(rpatch.norm(), sigma, epsPatches, aPatches, rstarPatches);
                    if (rpatch.norm() == 0) {
                        patchForce = vec3<double> (0, 0, 0);
                    }
                    else {
                        patchForce = patchesForceNorm*rpatch/rpatch.norm();
                    }
                    // Calculate force and torque acting on particle 1 and add values to previous forces and torques
                    force1 += patchForce;
                    torque1 += 0.5*sigma * patchNormal1.cross(patchForce);

                    // Calculate force and torque acting on particle 2 and add values to previous forces and torques
                    force2 += -1.0*patchForce;
                    torque2 += 0.5*sigma * patchNormal2.cross(-1.0*patchForce);
                }
            }
        }
        return {force + force1, torque1, -1.0*force + force2, torque2};
    }


}