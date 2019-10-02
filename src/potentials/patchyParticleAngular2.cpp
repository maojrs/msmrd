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
            vec3<double> derivativeAngularPotential = calculateQuaternionTorque(part1, part2);
            derivativeAngularPotential = angularStrength * derivativeAngularPotential;
            torque1 += derivativeAngularPotential;
            torque2 -= derivativeAngularPotential;
        }

        return {force + force1, torque1, -1.0 * force + force2, torque2};
    }

    /* Given two quaternions/orientations, returns a rotation in the axis-angle representation, representing
     * the torque that should be applied. */
    vec3<double> patchyParticleAngular2::calculateQuaternionTorque(particle part1, particle part2) {
        // Calculate relative orientation
        //quaternion<double> relativeOrientation = part2.orientation * part1.orientation.conj();

        // Calculate positions of patches
        vec3<double> part1Patch1 = part1.position +
                0.5*sigma * msmrdtools::rotateVec(patchesCoordinates[0], part1.orientation);
        vec3<double> part1Patch2 = part1.position +
                0.5*sigma * msmrdtools::rotateVec(patchesCoordinates[1], part1.orientation);
        vec3<double> part2Patch1 = part2.position +
                0.5*sigma * msmrdtools::rotateVec(patchesCoordinates[0], part2.orientation);
        vec3<double> part2Patch2 = part2.position +
                0.5*sigma * msmrdtools::rotateVec(patchesCoordinates[1], part2.orientation);
        // Calculate the 4 relevant distances between patches
        std::vector<double> patchesDistances(4);
        patchesDistances[0] = (part1Patch1 - part2Patch1).norm();
        patchesDistances[1] = (part1Patch1 - part2Patch2).norm();
        patchesDistances[2] = (part1Patch2 - part2Patch1).norm();
        patchesDistances[3] = (part1Patch2 - part2Patch2).norm();
        // Find minimum distance and match corresponding expected orientation (see setMetastableRegions())
        auto minIndex = std::min_element(patchesDistances.begin(), patchesDistances.end()) - patchesDistances.begin();

        // Calculate rotation of orthonormal basis (i,j,k)
        vec3<double> iUnitary = vec3<double>(1.0, 0.0, 0.0);
        vec3<double> jUnitary = vec3<double>(0.0, 1.0, 0.0);
        vec3<double> kUnitary = vec3<double>(0.0, 0.0, 1.0);
        // Obtain orthonormal basis fixed to particle 1
        vec3<double> iPart1 = msmrdtools::rotateVec(iUnitary, part1.orientation);
        vec3<double> jPart1 = msmrdtools::rotateVec(jUnitary, part1.orientation);
        vec3<double> kPart1 = msmrdtools::rotateVec(kUnitary, part1.orientation);
        // Obtain current orthonormal basis fixed to particle 2
        vec3<double> iPart2 = msmrdtools::rotateVec(iUnitary, part2.orientation);
        vec3<double> jPart2 = msmrdtools::rotateVec(jUnitary, part2.orientation);
        vec3<double> kPart2 = msmrdtools::rotateVec(kUnitary, part2.orientation);
        /* Obtain expected orthonormal basis of particle2, corresponds to an orientation
         * minima, i.e. the orientation we wish the system to relax to.*/
        vec3<double> iPart2Expected = msmrdtools::rotateVec(iPart1, quatRotations[minIndex]);
        vec3<double> jPart2Expected = msmrdtools::rotateVec(jPart1, quatRotations[minIndex]);
        vec3<double> kPart2Expected = msmrdtools::rotateVec(kPart1, quatRotations[minIndex]);

        /* The goal is to produce a torque that always brings the orthonormal basis of
         * particle 2, closer to the expected one. This can be simply done with the cross product
         * (equivalent to a potential of -cos(theta_i) or a more complicated one -1.0 *like -(cos(theta_i) + 1)^8,
         * where theta is the angle between the orthonormal basis of particle 2 and its expected one.)*/
        vec3<double> torque = vec3<double>(0.0, 0.0, 0.0);
        // Derivative of potential -(cos(theta) + 1)^8
        double cos1sqr = (iPart2 * iPart2Expected + 1) * (iPart2 * iPart2Expected + 1);
        double cos2sqr = (jPart2 * jPart2Expected + 1) * (jPart2 * jPart2Expected + 1);
        double cos3sqr = (kPart2 * kPart2Expected + 1) * (kPart2 * kPart2Expected + 1);
        torque += 8 * cos1sqr * cos1sqr * cos1sqr * (iPart2 * iPart2Expected + 1) * iPart2.cross(iPart2Expected);
        torque += 8 * cos2sqr * cos2sqr * cos2sqr * (jPart2 * jPart2Expected + 1) * jPart2.cross(jPart2Expected);
        torque += 8 * cos3sqr * cos3sqr * cos3sqr * (kPart2 * kPart2Expected + 1) * kPart2.cross(kPart2Expected);
//        torque += iPart2.cross(iPart2Expected);
//        torque += jPart2.cross(jPart2Expected);
//        torque += kPart2.cross(kPart2Expected);

        /* The minus one is to make it consistent with the sign choice in forceTorque function. */
        return -1.0 * torque;
    }


    /* Sets bound states (metastable regions) of this patchy dimer implementation. Same to setBoundStates
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