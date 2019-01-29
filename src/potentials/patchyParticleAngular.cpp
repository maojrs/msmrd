//
// Created by maojrs on 1/14/19.
//

#include <utility>
#include "potentials/patchyParticleAngular.hpp"


namespace msmrd {
    /*
     * Constructors inherited from patchyParticle. Added additional constructor to incorporate
     * angular potential strength
     */
    patchyParticleAngular::patchyParticleAngular(double sigma, double strength, double angularStrength,
                                                 std::vector<std::vector<double>> patchesCoordinates) :
    angularStrength(angularStrength), patchyParticle(sigma, strength, patchesCoordinates) {};

    patchyParticleAngular::patchyParticleAngular(double sigma, double strength, double angularStrength,
                                                 std::vector<vec3<double>> patchesCoordinates) :
            angularStrength(angularStrength), patchyParticle(sigma, strength, patchesCoordinates) {};


    // Evaluates potential at given positions and orientations of two particles
    double patchyParticleAngular::evaluate(const particle &part1, const particle &part2) {
        double repulsivePotential;
        double attractivePotential;
        double patchesPotential = 0.0;
        double angularPotential;
        vec3<double> patch1;
        vec3<double> patch2;
        vec3<double> patchNormal1;
        vec3<double> patchNormal2;
        vec3<double> rpatch;

        vec3<double> pos1 = part1.position;
        vec3<double> pos2 = part2.position;
        quaternion<double> theta1 = part1.orientation;
        quaternion<double> theta2 = part2.orientation;
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
                    patchesPotential += patchPotentialScaling*quadraticPotential(rpatch.norm(), sigma, epsPatches, aPatches, rstarPatches);
                }
            }
        }

        /* Explicit angular dependence based on first two patches of the two particles (not most efficient approach
         * but efficiency is not a problem in this example) */
        if (rvec.norm() <= 2*sigma) {
            // Calculate all normal vectors to first two patches for both particles
            vec3<double> part1PatchNormal1 = msmrdtools::rotateVec(patchesCoordinates[0], theta1);
            vec3<double> part1PatchNormal2 = msmrdtools::rotateVec(patchesCoordinates[1], theta1);
            vec3<double> part2PatchNormal1 = msmrdtools::rotateVec(patchesCoordinates[0], theta2);
            vec3<double> part2PatchNormal2 = msmrdtools::rotateVec(patchesCoordinates[1], theta2);

            // Calculate unitary vectors describing planes where particle center and first two patches are
            vec3<double> plane1 = part1PatchNormal1.cross(part1PatchNormal2);
            vec3<double> plane2 = part2PatchNormal1.cross(part2PatchNormal2);
            plane1 = plane1/plane1.norm();
            plane2 = plane2/plane2.norm();

            // Additional rotational potential depending on how well aligned the planes are. (-cos(theta)^8)
            double cosSquared = (plane1*plane2)*(plane1*plane2);
            angularPotential = -0.5*sigma*angularStrength*cosSquared*cosSquared*cosSquared*cosSquared;
        }


        return repulsivePotential + attractivePotential + patchesPotential + angularPotential;
    }

    /* Calculate and return (force1, torque1, force2, torque2), which correspond to the force and torque
     * acting on particle1 and the force and torque acting on particle2, respectively. */
    std::array<vec3<double>, 4> patchyParticleAngular::forceTorque(const particle &part1, const particle &part2) {
        vec3<double> force;
        vec3<double> force1 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> force2 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> torque1 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> torque2 = vec3<double> (0.0, 0.0, 0.0);

        vec3<double> pos1 = part1.position;
        vec3<double> pos2 = part2.position;
        quaternion<double> theta1 = part1.orientation;
        quaternion<double> theta2 = part2.orientation;
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
                        patchForce = patchesForceNorm*rpatch/rpatch.norm();
                    }
                    // Calculate force and torque acting on particle 1 and add values to previous forces and torques
                    force1 += patchPotentialScaling * patchForce;
                    torque1 += 0.5*sigma * patchPotentialScaling * patchNormal1.cross(patchForce);

                    // Calculate force and torque acting on particle 2 and add values to previous forces and torques
                    force2 += -1.0 * patchPotentialScaling * patchForce;
                    torque2 += 0.5*sigma * patchPotentialScaling * patchNormal2.cross(-1.0*patchForce);
                }
            }
        }

        /* Explicit angular dependence based on first two patches of the two particles (not most efficient approach
         * but efficiency is not a problem in this example) */
        if (rvec.norm() <= 2*sigma) {
            // Calculate all normal vectors to first two patches for both particles
            vec3<double> part1PatchNormal1 = msmrdtools::rotateVec(patchesCoordinates[0], theta1);
            vec3<double> part1PatchNormal2 = msmrdtools::rotateVec(patchesCoordinates[1], theta1);
            vec3<double> part2PatchNormal1 = msmrdtools::rotateVec(patchesCoordinates[0], theta2);
            vec3<double> part2PatchNormal2 = msmrdtools::rotateVec(patchesCoordinates[1], theta2);

            // Calculate unitary vectors describing planes where particle center and first two patches are
            vec3<double> plane1 = part1PatchNormal1.cross(part1PatchNormal2);
            vec3<double> plane2 = part2PatchNormal1.cross(part2PatchNormal2);
            plane1 = plane1/plane1.norm();
            plane2 = plane2/plane2.norm();

            // Additional torques depending on how well aligned the planes are. (Derivative of -cos(theta)^8)
            double cosSquared = (plane1*plane2)*(plane1*plane2);
            double cosSeventh = cosSquared*cosSquared*cosSquared*(plane1*plane2);
            vec3<double> derivativeAngluarPotential = 0.5*angularStrength*sigma*8*cosSeventh*plane1.cross(plane2);
            torque1 += derivativeAngluarPotential;
            torque2 -= derivativeAngluarPotential;
        }

        return {force + force1, torque1, -1.0*force + force2, torque2};
    }

}