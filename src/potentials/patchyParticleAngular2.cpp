//
// Created by maojrs on 9/4/19.
//

#include "potentials/patchyParticleAngular2.hpp"
#include "tools.hpp"


namespace msmrd {
    /*
     * Constructors inherited from patchyParticleAngular.
     */
    patchyParticleAngular2::patchyParticleAngular2(double sigma, double strength, double angularStrength,
                           std::vector<std::vector<double>> patchesCoordinates) :
            patchyParticleAngular(sigma, strength, angularStrength, patchesCoordinates) {
        setMetastableRegions();
    };

    patchyParticleAngular2::patchyParticleAngular2(double sigma, double strength, double angularStrength,
                                                   std::vector<vec3<double>> patchesCoordinates) :
            patchyParticleAngular(sigma, strength, angularStrength, patchesCoordinates) {
        setMetastableRegions();
    };



    // Evaluates potential at given positions and orientations of two particles
    double patchyParticleAngular2::evaluate(const particle &part1, const particle &part2) {
          // Not available for force calculated using quaternions (at least not yet).
        return 0;
    }

    /* Calculate and return (force1, torque1, force2, torque2), which correspond to the force and torque
     * acting on particle1 and the force and torque acting on particle2, respectively. */
    std::array<vec3<double>, 4> patchyParticleAngular2::forceTorque(const particle &part1, const particle &part2) {
        vec3<double> force;
        vec3<double> force1 = vec3<double>(0.0, 0.0, 0.0);
        vec3<double> force2 = vec3<double>(0.0, 0.0, 0.0);
        vec3<double> torque1 = vec3<double>(0.0, 0.0, 0.0);
        vec3<double> torque2 = vec3<double>(0.0, 0.0, 0.0);

        vec3<double> pos1 = part1.position;
        vec3<double> pos2 = part2.position;
        quaternion<double> theta1 = part1.orientation;
        quaternion<double> theta2 = part2.orientation;
        std::array<vec3<double>, 2> relPos = relativePositionComplete(pos1, pos2);
        vec3<double> pos1virtual = relPos[0]; // virtual pos1 if periodic boundary; otherwise pos1.
        vec3<double> rvec = relPos[1]; //pos2 - pos1;

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
        attractiveForceNorm = derivativeQuadraticPotential(rvec.norm(), sigma, epsAttractive, aAttractive,
                                                           rstarAttractive);
        force = (repulsiveForceNorm + attractiveForceNorm) * rvec / rvec.norm();

        // Calculate forces and torque due to patches interaction
        if (rvec.norm() <= 2 * sigma) {
            // Loop over all patches of particle 1
            for (int i = 0; i < patchesCoordinates.size(); i++) {
                patchNormal1 = msmrdtools::rotateVec(patchesCoordinates[i], theta1);
                patchNormal1 = patchNormal1 / patchNormal1.norm();
                patch1 = pos1virtual + 0.5 * sigma * patchNormal1;
                // Loop over all patches of particle 2
                for (int j = 0; j < patchesCoordinates.size(); j++) {
                    patchNormal2 = msmrdtools::rotateVec(patchesCoordinates[j], theta2);
                    patchNormal2 = patchNormal2 / patchNormal2.norm();
                    patch2 = pos2 + 0.5 * sigma * patchNormal2;
                    rpatch = patch2 - patch1;
                    // Calculate force vector between patches , correct sign of force given by rpatch/rpatch.norm().
                    patchesForceNorm = derivativeQuadraticPotential(rpatch.norm(), sigma, epsPatches, aPatches,
                                                                    rstarPatches);
                    if (rpatch.norm() == 0) {
                        patchForce = vec3<double>(0, 0, 0);
                    } else {
                        patchForce = patchesForceNorm * rpatch / rpatch.norm();
                    }
                    // Calculate force and torque acting on particle 1 and add values to previous forces and torques
                    force1 += patchPotentialScaling * patchForce;
                    torque1 += 0.5 * sigma * patchPotentialScaling * patchNormal1.cross(patchForce);

                    // Calculate force and torque acting on particle 2 and add values to previous forces and torques
                    force2 += -1.0 * patchPotentialScaling * patchForce;
                    torque2 += 0.5 * sigma * patchPotentialScaling * patchNormal2.cross(-1.0 * patchForce);
                }
            }
        }

        /* Explicit angular dependence using quaternions. */
        if (rvec.norm() <= 1.5 * sigma) {
            vec3<double> derivativeAngularPotential = calculateTorque(part1, part2);
            derivativeAngularPotential = angularStrength * derivativeAngularPotential;
            torque1 += derivativeAngularPotential;
            torque2 -= derivativeAngularPotential;
        }

        return {force + force1, torque1, -1.0 * force + force2, torque2};
    }


    /* Given two quaternions/orientations, returns a torque that should be applied. */
    vec3<double> patchyParticleAngular2::calculateTorque(particle part1, particle part2) {
        // Calculate relative orientation
        //quaternion<double> relativeOrientation = part2.orientation * part1.orientation.conj();

        // Calculate all normal vectors to first two patches for both particles
        vec3<double> part1PatchNormal1 = msmrdtools::rotateVec(patchesCoordinates[0], part1.orientation);
        vec3<double> part1PatchNormal2 = msmrdtools::rotateVec(patchesCoordinates[1], part1.orientation);
        vec3<double> part2PatchNormal1 = msmrdtools::rotateVec(patchesCoordinates[0], part2.orientation);
        vec3<double> part2PatchNormal2 = msmrdtools::rotateVec(patchesCoordinates[1], part2.orientation);

        // Calculate positions of patches
        vec3<double> part1Patch1 = part1.position + 0.5*sigma * part1PatchNormal1;
        vec3<double> part1Patch2 = part1.position + 0.5*sigma * part1PatchNormal2;
        vec3<double> part2Patch1 = part2.position + 0.5*sigma * part2PatchNormal1;
        vec3<double> part2Patch2 = part2.position + 0.5*sigma * part2PatchNormal2;
        // Calculate the 4 relevant distances between patches
        std::vector<double> patchesDistances(4);
        patchesDistances[0] = (part1Patch1 - part2Patch1).norm();
        patchesDistances[1] = (part1Patch1 - part2Patch2).norm();
        patchesDistances[2] = (part1Patch2 - part2Patch1).norm();
        patchesDistances[3] = (part1Patch2 - part2Patch2).norm();
        // Find minimum distance and match corresponding expected orientation (see setMetastableRegions())
        auto minIndex = static_cast<int> (std::min_element(patchesDistances.begin(), patchesDistances.end()) -
                                          patchesDistances.begin());

        // Calculate unitary vectors describing planes where particle center and first two patches are
        vec3<double> plane1;
        vec3<double> plane2;
        if (minIndex == 0 or minIndex == 3) {
            // Bound thorugh pathches (1,1) or (2,2)
            plane1 = part1PatchNormal1.cross(part1PatchNormal2);
            plane2 = part2PatchNormal1.cross(part2PatchNormal2);
        } else {
            // Bound thorugh pathches (1,2) or (2,1)
            plane1 = part1PatchNormal1.cross(part1PatchNormal2);
            plane2 = part2PatchNormal2.cross(part2PatchNormal1);
        }

        plane1 = plane1/plane1.norm();
        plane2 = plane2/plane2.norm();

        // Additional torques depending on how well aligned the planes are.
        vec3<double> torque;
        // It one minima, use potential of -(1/250)*(cos(theta) + 1)^8
        double cosSquared = (plane1 * plane2 + 1) * (plane1 * plane2 + 1);
        double cosSeventh = cosSquared * cosSquared * cosSquared * (plane1 * plane2 + 1);
        torque = (8.0/250.0) * cosSeventh * plane1.cross(plane2);
        return -1.0*torque;
    }


    /* DEPRECATED, no longer needed by calculateTorque, left for reference of metastable regions or future uses.
     * Sets bound states (metastable regions) of this patchy dimer implementation. Same to setBoundStates
     * in patchyDimer2 in trajectory/discrete/patchyDimer. */
    void patchyParticleAngular2::setMetastableRegions() {
        double angleDiff = 3 * M_PI / 5; // angle difference to form a pentamer
        /* Define relative position vectors from particle 1 at the origin. These two patches
         * point in the same direction as the two patches in the dimer. */
        refRelativePositions[0] = {std::cos(angleDiff / 2.0), std::sin(angleDiff / 2.0), 0};
        refRelativePositions[1] = {std::cos(angleDiff / 2.0), std::sin(-angleDiff / 2.0), 0};
        vec3<double> relPos1orthogonal = {-1.0 * std::sin(angleDiff / 2.0), std::cos(angleDiff / 2.0), 0.0};
        vec3<double> relPos2orthogonal = {std::sin(angleDiff / 2.0), std::cos(angleDiff / 2.0), 0.0};
        /* Relative rotations (from particle 1) of particle 2 that yield the 4 bound states
         * in the axis-angle representation. (One needs to make drawing to understand)*/
        std::array<vec3<double>, 4> rotations;
        rotations[0] = M_PI * relPos1orthogonal; // part1Patch1 with part2Patch1 (see calculateQuaternionTorque)
        rotations[1] = {0.0, 0.0, -2 * M_PI / 5.0}; // part1Patch1 with part2Patch2
        // --first 2 rotations correspond to binding on top patch of particle 1, next 2 rotations to bottom patch
        rotations[2] = {0.0, 0.0, 2 * M_PI / 5.0}; // part1Patch2 with part2Patch1
        rotations[3] = M_PI * relPos2orthogonal; // part1Patch2 with part2Patch2
        for (int i = 0; i < 4; i++) {
            quatRotations[i] = msmrdtools::axisangle2quaternion(rotations[i]);
        }
    }

}