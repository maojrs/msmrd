//
// Created by maojrs on 10/23/19.
//

#include "potentials/patchyProteinAngular.hpp"

namespace msmrd{

    // Evaluates potential at given positions and orientations of two particles
    double patchyProteinAngular::evaluate(const particleMS &part1, const particleMS &part2) {
        // Decalre variables used in loop
        double patchesPotential = 0.0;
        vec3<double> patch1;
        vec3<double> patch2;
        vec3<double> patchNormal1;
        vec3<double> patchNormal2;
        vec3<double> rpatch;

        // Calculates relative position
        auto relPos = relativePositionComplete(part1.position, part2.position);
        vec3<double> pos1virtual = relPos[0]; // virtual part1.position if periodic boundary; otherwise part1.position.
        vec3<double> rvec = relPos[1]; // part2.position - part1.position;

        // Calculate isotropic potential
        auto repulsivePotential = quadraticPotential(rvec.norm(), sigma, epsRepulsive, aRepulsive, rstarRepulsive);
        auto attractivePotential = quadraticPotential(rvec.norm(), sigma, epsAttractive, aAttractive, rstarAttractive);

        /* Assign patch pattern depending on particle type (note only two types of particles are supported here) */
        auto patchesCoords1 = assignPatches(part1.type);
        auto patchesCoords2 = assignPatches(part2.type);

        // WORKING HERE

        // Check for patches interactions, only if particle 2 is in state 0
        if ((rvec.norm() <= 2*sigma) and part2.state == 1) {
            // Loop over all patches
            for (int i = 0; i < patchesCoords1.size(); i++) {
                patchNormal1 = msmrdtools::rotateVec(patchesCoords1[i], part1.orientation);
                patch1 = pos1virtual + 0.5*sigma*patchNormal1;
                for (int j = 0; j < patchesCoords2.size(); j++) {
                    patchNormal2 = msmrdtools::rotateVec(patchesCoords2[j], part2.orientation);
                    patch2 = part2.position + 0.5*sigma*patchNormal2;
                    rpatch = patch2 - patch1; // Scale unit distance of patches by sigma
                    // Assumes the first patch from type1 has a different type of interaction,
                    if ( (i == 0 && part1.type == 0) || (j == 0 && part2.type == 0) ) {
                        patchesPotential += quadraticPotential(rpatch.norm(), sigma, epsPatches[1],
                                                               aPatches[1], rstarPatches[1]);
                    }
                        // while all the other patches combinations have another type of interaction.
                    else{
                        patchesPotential += quadraticPotential(rpatch.norm(), sigma, epsPatches[0],
                                                               aPatches[0], rstarPatches[0]);
                    }
                }
            }
        }

        return repulsivePotential + attractivePotential + patchesPotential;
    }

    /* Calculate and return (force1, torque1, force2, torque2), which correspond to the force and torque
     * acting on particle1 and the force and torque acting on particle2, respectively. */
    std::array<vec3<double>, 4> patchyProteinAngular::forceTorque(const particleMS &part1, const particleMS &part2)  {

        vec3<double> force1 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> force2 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> torque1 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> torque2 = vec3<double> (0.0, 0.0, 0.0);

        // Calculate relative position
        std::array<vec3<double>, 2> relPos = relativePositionComplete(part1.position, part2.position);
        vec3<double> pos1virtual = relPos[0]; // virtual part1.position if periodic boundary; otherwise part1.position.
        vec3<double> rvec = relPos[1]; //part2.position - part1.position;

        // auxiliary variables to calculate force and torque
        double patchesForceNorm;
        vec3<double> patchForce;
        vec3<double> patch1;
        vec3<double> patch2;
        vec3<double> patchNormal1;
        vec3<double> patchNormal2;
        vec3<double> rpatch;

        /* Calculate and add forces due to repulsive and attractive isotropic potentials.
         *  Note correct sign/direction of force given by rvec/rvec.norm*() */
        auto repulsiveForceNorm = derivativeQuadraticPotential(rvec.norm(), sigma, epsRepulsive,
                                                               aRepulsive, rstarRepulsive);
        auto attractiveForceNorm = derivativeQuadraticPotential(rvec.norm(), sigma, epsAttractive,
                                                                aAttractive, rstarAttractive);
        auto force = (repulsiveForceNorm + attractiveForceNorm)*rvec/rvec.norm();

        /* Assign patch pattern depending on particle type (note only two types of particles are supported here) */
        auto patchesCoords1 = assignPatches(part1.type);
        auto patchesCoords2 = assignPatches(part2.type);

        // Calculate forces and torque due to patches interaction
        if (rvec.norm() <= 2*sigma) {
            // Loop over all patches of particle 1
            for (int i = 0; i < patchesCoords1.size(); i++) {
                patchNormal1 = msmrdtools::rotateVec(patchesCoords1[i], part1.orientation);
                patchNormal1 = patchNormal1/patchNormal1.norm();
                patch1 = pos1virtual + 0.5*sigma*patchNormal1;
                // Loop over all patches of particle 2
                for (int j = 0; j < patchesCoords2.size(); j++) {
                    patchNormal2 = msmrdtools::rotateVec(patchesCoords2[j], part2.orientation);
                    patchNormal2 = patchNormal2/patchNormal2.norm();
                    patch2 = part2.position + 0.5*sigma*patchNormal2;
                    // Calculate distance between the two patches
                    rpatch = patch2 - patch1;
                    /* Calculate force vector between patches , correct sign of force given by rpatch/rpatch.norm().
                     * It also assumes the first patch from type = 0 has a different type of interaction. */
                    if ( (i == 0 && part1.type == 0) || (j == 0 && part2.type == 0) ) {
                        patchesForceNorm = derivativeQuadraticPotential(rpatch.norm(), sigma, epsPatches[1],
                                                                        aPatches[1], rstarPatches[1]);
                    }
                    else {
                        patchesForceNorm = derivativeQuadraticPotential(rpatch.norm(), sigma, epsPatches[0],
                                                                        aPatches[0], rstarPatches[0]);
                    }
                    // Determine force vector avoiding division by zero
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


    /* Given two quaternions/orientations, returns planes(unit vectors) to be aligned by torque. These
     * may vary depending on the physical arrangement of your molecules and how the patches are ordered and
     * selected. Currently set up for particle 1 with 6 binding patches and particle 2 only with one.*/
    std::tuple<vec3<double>, vec3<double>> patchyProteinAngular::calculatePlanes(const particleMS &part1,
                                                                 const particleMS &part2,
                                                                 const std::vector<vec3<double>> patches1,
                                                                 const std::vector<vec3<double>> patches2) {

        std::vector<vec3<double>> part1PatchNormals;
        std::vector<vec3<double>> part2PatchNormals;
        vec3<double> patchRef = vec3<double> (0.0, 0.0, 1.0);
        for (auto &patch : patches1) {
            part1PatchNormals.push_back(msmrdtools::rotateVec(patch, part1.orientation));
        }
        for (auto &patch : patches2) {
            part2PatchNormals.push_back(msmrdtools::rotateVec(patch, part2.orientation));
        }

        // Calculate positions of patches
        std::vector<vec3<double>> part1PatchPositions;
        std::vector<vec3<double>> part2PatchPositions;
        for(auto &patchNormal : part1PatchNormals) {
            part1PatchPositions.push_back(part1.position + 0.5*sigma * patchNormal);
        }
        for(auto &patchNormal : part2PatchNormals) {
            part2PatchPositions.push_back(part2.position + 0.5*sigma * patchNormal);
        }


        /* Calculate all the relevant distances between patches (initially coded for 6 patches in particle1 and
         * 1 patch in particle2) */
        std::vector<double> patchesDistances;
        for (auto part1Patch : part1PatchPositions) {
            for (auto part2Patch : part2PatchPositions) {
                patchesDistances.push_back((part1Patch - part2Patch).norm());
            }
        }

        /* Find minimum distance and match corresponding expected orientation (0, 1, 2 ... 7 in original implementation)
         * depending on the patches that are closer to each other) */
        auto minIndex = static_cast<int> (std::min_element(patchesDistances.begin(), patchesDistances.end()) -
                                          patchesDistances.begin());

        // Calculate unitary vectors describing planes where particle center and first two patches are
        vec3<double> plane1;
        vec3<double> plane2;

        // Planes to align  selection dependent on implementation
        int secondIndex = (minIndex + 1) % 8; // 8=part1PatchNormals.size()
        patchRef = msmrdtools::rotateVec(patchRef, part2.orientation);
        plane1 = part1PatchNormals[minIndex].cross(part1PatchNormals[secondIndex]); // should always be perpendicular.
        plane2 = part2PatchNormals[0].cross(patchRef);

        return std::make_tuple(plane1, plane2);
    }


}

