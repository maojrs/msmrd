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
    angularStrength(angularStrength), patchyParticle(sigma, strength, patchesCoordinates) {
        patchPotentialScaling = 0.5; // 0.3 original value
    };

    patchyParticleAngular::patchyParticleAngular(double sigma, double strength, double angularStrength,
                                                 std::vector<vec3<double>> patchesCoordinates) :
            angularStrength(angularStrength), patchyParticle(sigma, strength, patchesCoordinates) {
        patchPotentialScaling = 0.5; // 0.3 original value
    };


    // Evaluates potential at given positions and orientations of two particles
    double patchyParticleAngular::evaluate(particle &part1, particle &part2) {

        // Get part of potential that is the same as for normal patchy particle from parent function.
        double patchyParticlePotential = patchyParticle::evaluate(part1, part2);

        std::array<vec3<double>, 2> relPos = relativePositionComplete(part1.position, part2.position);
        vec3<double> rvec = relPos[1]; //part2.position - part1.position;

        /* Explicit angular dependence based on first two patches of the two particles (not most efficient approach
         * but efficiency is not a problem in this example) */
        double angularPotential = 0.0;
        if (rvec.norm() <= 2*sigma and patchesActive) {
            // Calculate all normal vectors to first two patches for both particles
            vec3<double> part1PatchNormal1 = msmrdtools::rotateVec(patchesCoordinates[0], part1.orientation);
            vec3<double> part1PatchNormal2 = msmrdtools::rotateVec(patchesCoordinates[1], part1.orientation);
            vec3<double> part2PatchNormal1 = msmrdtools::rotateVec(patchesCoordinates[0], part2.orientation);
            vec3<double> part2PatchNormal2 = msmrdtools::rotateVec(patchesCoordinates[1], part2.orientation);

            // Calculate unitary vectors describing planes where particle center and first two patches are
            vec3<double> plane1 = part1PatchNormal1.cross(part1PatchNormal2);
            vec3<double> plane2 = part2PatchNormal1.cross(part2PatchNormal2);
            plane1 = plane1/plane1.norm();
            plane2 = plane2/plane2.norm();

            /* Additional rotational potential depending on how well aligned
             * the planes are, (-cos(theta)^8) it has two minima at 0 and pi. */
            double cosSquared = (plane1 * plane2) * (plane1 * plane2);
            angularPotential = -0.5 * sigma * angularStrength * cosSquared * cosSquared * cosSquared * cosSquared;
        }


        return patchyParticlePotential + angularPotential;
    }

    /* Calculate and return (force1, torque1, force2, torque2), which correspond to the force and torque
     * acting on particle1 and the force and torque acting on particle2, respectively. */
    std::array<vec3<double>, 4> patchyParticleAngular::forceTorque(particle &part1, particle &part2) {

        // Get part of force and torque that is the same as for normal patchy particle from parent function.
        auto patchyParticleForceTorque = patchyParticle::forceTorque(part1, part2);

        vec3<double> force1 = patchyParticleForceTorque[0];
        vec3<double> torque1 = patchyParticleForceTorque[1];
        vec3<double> force2 = patchyParticleForceTorque[2];
        vec3<double> torque2 = patchyParticleForceTorque[3];

        std::array<vec3<double>, 2> relPos = relativePositionComplete(part1.position, part2.position);
        vec3<double> pos1virtual = relPos[0]; // virtual pos1 if periodic boundary; otherwise pos1.
        vec3<double> rvec = relPos[1]; //pos2 - pos1;

        /* Explicit angular dependence based on first two patches of the two particles (not most efficient approach
         * but efficiency is not a problem in this example) */
        if (rvec.norm() <= 2*sigma and patchesActive) {
            // Calculate all normal vectors to first two patches for both particles
            vec3<double> part1PatchNormal1 = msmrdtools::rotateVec(patchesCoordinates[0], part1.orientation);
            vec3<double> part1PatchNormal2 = msmrdtools::rotateVec(patchesCoordinates[1], part1.orientation);
            vec3<double> part2PatchNormal1 = msmrdtools::rotateVec(patchesCoordinates[0], part2.orientation);
            vec3<double> part2PatchNormal2 = msmrdtools::rotateVec(patchesCoordinates[1], part2.orientation);

            // Calculate unitary vectors describing planes where particle center and first two patches are
            vec3<double> plane1 = part1PatchNormal1.cross(part1PatchNormal2);
            vec3<double> plane2 = part2PatchNormal1.cross(part2PatchNormal2);
            plane1 = plane1/plane1.norm();
            plane2 = plane2/plane2.norm();

            // Torque calculation depending on how well aligned the planes are.
            vec3<double> derivativeAngluarPotential;
            // It two minimas, use potential of -cos(theta)^8
            double cosSquared = (plane1 * plane2) * (plane1 * plane2);
            double cosSeventh = cosSquared * cosSquared * cosSquared * (plane1 * plane2);
            derivativeAngluarPotential = 0.5 * angularStrength * sigma * 8 * cosSeventh * plane1.cross(plane2);
            torque1 += derivativeAngluarPotential; // Plus sign since plane1 x plane2 defined torque in particle 1
            torque2 -= derivativeAngluarPotential;
        }

        return {force1, torque1, force2, torque2};
    }



    /*
     * Begins implementations of second version of patchy particle angular potential, ie patchyParticleAngular2
     */

    // Evaluates potential at given positions and orientations of two particles
    double patchyParticleAngular2::evaluate(particle &part1, particle &part2) {

        // Get part of potential that is the same as for normal patchy particle from parent function.
        double patchyParticlePotential = patchyParticle::evaluate(part1, part2);

        std::array<vec3<double>, 2> relPos = relativePositionComplete(part1.position, part2.position);
        vec3<double> rvec = relPos[1]; //part2.position - part1.position;

        /* Explicit angular dependence based on first two patches of the two particles (not most efficient approach
         * but efficiency is not a problem in this example). Alternative version as patchyParticleAngular */
        double angularPotential = 0.0;
        if (rvec.norm() <= 2.0 * sigma and patchesActive) {
            /* Get planes needed to be aligned by torque, based on use potential of -[(cos(theta) + 1)/2]^8
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
    std::array<vec3<double>, 4> patchyParticleAngular2::forceTorque(particle &part1, particle &part2) {

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
        vec3<double> derivativeAngluarPotential = vec3<double> (0.0, 0.0, 0.0);
        if (rvec.norm() <= 2.0 * sigma and patchesActive) {
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
    std::tuple<vec3<double>, vec3<double>> patchyParticleAngular2::calculatePlanes(const particle &part1,
                                                                                   const particle &part2) {

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

}
