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
    double patchyParticleAngular::evaluate(const particle &part1, const particle &part2) {

        // Get part of potential that is the same as for normal patchy particle from parent function.
        double patchyParticlePotential = patchyParticle::evaluate(part1, part2);

        std::array<vec3<double>, 2> relPos = relativePositionComplete(part1.position, part2.position);
        vec3<double> rvec = relPos[1]; //part2.position - part1.position;

        /* Explicit angular dependence based on first two patches of the two particles (not most efficient approach
         * but efficiency is not a problem in this example) */
        double angularPotential = 0.0;
        if (rvec.norm() <= 2*sigma) {
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
    std::array<vec3<double>, 4> patchyParticleAngular::forceTorque(const particle &part1, const particle &part2) {

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
        if (rvec.norm() <= 2*sigma) {
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

}