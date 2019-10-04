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

        // Get part of potential that is the same as for normal patchy particle from parent function.
        double patchyParticlePotential = patchyParticle::evaluate(part1, part2);

        std::array<vec3<double>, 2> relPos = relativePositionComplete(part1.position, part2.position);
        vec3<double> rvec = relPos[1]; //part2.position - part1.position;

        /* Explicit angular dependence based on first two patches of the two particles (not most efficient approach
         * but efficiency is not a problem in this example). Alternative version as patchyParticleAngular */
        double angularPotential = 0.0;
        if (rvec.norm() <= 2.0 * sigma) {
            /* Get planes needed to be aligned by torque, based on use potential of -(1/256)*(cos(theta) + 1)^8
             * with only one minima */
            vec3<double> plane1;
            vec3<double> plane2;
            std::tie(plane1, plane2) = calculatePlanes(part1, part2);
            double cosSquared = (plane1 * plane2 + 1) * (plane1 * plane2 + 1);
            angularPotential = - (1.0/256.0) * angularStrength * cosSquared * cosSquared * cosSquared * cosSquared;
        }

        return patchyParticlePotential + angularPotential;

    }

    /* Calculate and return (force1, torque1, force2, torque2), which correspond to the force and torque
     * acting on particle1 and the force and torque acting on particle2, respectively. */
    std::array<vec3<double>, 4> patchyParticleAngular2::forceTorque(const particle &part1, const particle &part2) {

        // Get part of force and torque that is the same as for normal patchy particle from parent function.
        auto patchyParticleForceTorque = patchyParticle::forceTorque(part1, part2);

        vec3<double> force1 = patchyParticleForceTorque[0];
        vec3<double> torque1 = patchyParticleForceTorque[1];
        vec3<double> force2 = patchyParticleForceTorque[2];
        vec3<double> torque2 = patchyParticleForceTorque[3];

        std::array<vec3<double>, 2> relPos = relativePositionComplete(part1.position, part2.position);
        vec3<double> pos1virtual = relPos[0]; // virtual pos1 if periodic boundary; otherwise pos1.
        vec3<double> rvec = relPos[1]; //pos2 - pos1;

        /* Explicit angular dependence. */
        vec3<double> derivativeAngluarPotential;
        if (rvec.norm() <= 2.0 * sigma) {
            // Get planes needed to be aligned by torque
            vec3<double> plane1;
            vec3<double> plane2;
            std::tie(plane1, plane2) = calculatePlanes(part1, part2);
            /* If one uses only one minima, use potential of -(1/256)*(cos(theta) + 1)^8, with theta the angle
             * between the unitary vector of each plane */
            double cosSquared = (plane1 * plane2 + 1) * (plane1 * plane2 + 1);
            double cosSeventh = cosSquared * cosSquared * cosSquared * (plane1 * plane2 + 1);
            derivativeAngluarPotential = 0.5 * angularStrength * sigma * (8.0/256.0) * cosSeventh * plane1.cross(plane2);
            torque1 += derivativeAngluarPotential; // Plus sign since plane1 x plane2 defined torque in particle 1
            torque2 -= derivativeAngluarPotential;
        }

        return {force1, torque1, force2, torque2};
    }


    /* Given two quaternions/orientations, returns planes(unit vectors) to be aligned by torque. These
     * may vary depending on the physical arrangement of your molecules. */
    std::tuple<vec3<double>, vec3<double>> patchyParticleAngular2::calculatePlanes(particle part1, particle part2) {

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

        // Calculate the 4 relevant distances between patches (1,1), (1,2), (2,1) and (2,2)
        std::vector<double> patchesDistances(4);
        patchesDistances[0] = (part1Patch1 - part2Patch1).norm();
        patchesDistances[1] = (part1Patch1 - part2Patch2).norm();
        patchesDistances[2] = (part1Patch2 - part2Patch1).norm();
        patchesDistances[3] = (part1Patch2 - part2Patch2).norm();

        /* Find minimum distance and match corresponding expected orientation (0, 1, 2 or 3,
         * depending on the patches that are closer to each other) */
        auto minIndex = static_cast<int> (std::min_element(patchesDistances.begin(), patchesDistances.end()) -
                                          patchesDistances.begin());

        // Calculate unitary vectors describing planes where particle center and first two patches are
        vec3<double> plane1;
        vec3<double> plane2;
        if (minIndex == 1 or minIndex == 2) {
            // Bound thorugh pathches (1,2) or (2,1)
            plane1 = part1PatchNormal1.cross(part1PatchNormal2);
            plane2 = part2PatchNormal1.cross(part2PatchNormal2);
        } else {
            // Bound thorugh pathches (1,1) or (2,2)
            plane1 = part1PatchNormal1.cross(part1PatchNormal2);
            plane2 = part2PatchNormal2.cross(part2PatchNormal1);
        }

        plane1 = plane1/plane1.norm();
        plane2 = plane2/plane2.norm();

        return std::make_tuple(plane1, plane2);
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