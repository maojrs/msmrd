//
// Created by maojrs on 10/23/19.
//

#include <utility>
#include "potentials/patchyProteinMarkovSwitch.hpp"

namespace msmrd{
    /*
     * Constructor inherited from patchyProtein. This additional constructors are to incorporate
     * angular potential strength
     */
    patchyProteinMarkovSwitch::patchyProteinMarkovSwitch(double sigma, double strength, double angularStrength,
                                 std::vector<vec3<double>> patchesCoordinatesA,
                                 std::vector<vec3<double>> patchesCoordinatesB)
            :  patchyProtein(sigma, strength, std::move(patchesCoordinatesA), std::move(patchesCoordinatesB)),
               angularStrength(angularStrength) {};

    patchyProteinMarkovSwitch::patchyProteinMarkovSwitch(double sigma, double strength, double angularStrength,
                                 std::vector<std::vector<double>> patchesCoordsA,
                                 std::vector<std::vector<double>> patchesCoordsB)
            :  patchyProtein(sigma, strength, patchesCoordinatesA, patchesCoordinatesB),
               angularStrength(angularStrength) {};

    // Evaluates potential at given positions and orientations of two particles
    double patchyProteinMarkovSwitch::evaluate(const particleMS &part1, const particleMS &part2) {
        // Declare variables used in loop
        double patchesPotential = 0.0;
        double angularPotential = 0.0;

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

        // Evaluate patches potential
        if (part2.state == 1) {
            // Use default patches auxiliary parent function without angular dependence
            patchesPotential = evaluatePatchesPotential(part1, part2, relPos, patchesCoords1, patchesCoords2);
            /* Get planes needed to be aligned by torque, based on use potential of -[(cos(theta) + 1)/2]^8
             * with only one minima. This adds the angular dependence based on planes calculated in calculatePlanes.
             * Implementation specific.*/
            vec3<double> plane1;
            vec3<double> plane2;
            std::tie(plane1, plane2) = calculatePlanes(part1, part2, patchesCoords1, patchesCoords2);
            double cosSquared = (plane1 * plane2 + 1) * (plane1 * plane2 + 1);
            angularPotential = -(1.0 / 256.0) * angularStrength * cosSquared * cosSquared * cosSquared * cosSquared;
        }

        return repulsivePotential + attractivePotential + patchesPotential + angularPotential;
    }

    /* Calculate and return (force1, torque1, force2, torque2), which correspond to the force and torque
     * acting on particle1 and the force and torque acting on particle2, respectively. */
    std::array<vec3<double>, 4> patchyProteinMarkovSwitch::forceTorque(particleMS &part1, particleMS &part2)  {

        vec3<double> force1 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> force2 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> torque1 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> torque2 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> derivativeAngluarPotential = vec3<double> (0.0, 0.0, 0.0);

        // Calculate relative position
        std::array<vec3<double>, 2> relPos = relativePositionComplete(part1.position, part2.position);
        vec3<double> pos1virtual = relPos[0]; // virtual part1.position if periodic boundary; otherwise part1.position.
        vec3<double> rvec = relPos[1]; //part2.position - part1.position;

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

        /* Enables/Disables MSM following potential hardcoded rules. In this case, if bounded or close to bounded
         * disable the MSM on particles type 1 and state 0. Not always necessary, but left commented just in case. */
        //enableDisableMSM(rvec, part1, part2);

        // Calculate forces and torque due to patches interaction, if close enough and if particle 2 is in state 1
        if (part2.state == 1) {
            // Calculate forces and torque due to patches interaction using auxiliary function
            auto forcTorqPatches = forceTorquePatches(part1, part2, relPos, patchesCoords1, patchesCoords2);
            force1 = forcTorqPatches[0];
            torque1 = forcTorqPatches[1];
            force2 = forcTorqPatches[2];
            torque2 = forcTorqPatches[3];

            /* Explicit angular dependence. */
            vec3<double> derivativeAngluarPotential;
            if (rvec.norm() <= 2.0 * sigma) {
                // Get planes needed to be aligned by torque
                vec3<double> plane1;
                vec3<double> plane2;
                std::tie(plane1, plane2) = calculatePlanes(part1, part2, patchesCoords1, patchesCoords2);
                /* If one uses only one minima, use potential of -(1/256)*(cos(theta) + 1)^8, with theta the angle
                 * between the unitary vector of each plane */
                double cosSquared = (plane1 * plane2 + 1) * (plane1 * plane2 + 1);
                double cosSeventh = cosSquared * cosSquared * cosSquared * (plane1 * plane2 + 1);
                derivativeAngluarPotential = 0.5 * angularStrength * sigma * (8.0/256.0) * cosSeventh * plane1.cross(plane2);
                torque1 += derivativeAngluarPotential; // Plus sign since plane1 x plane2 defined torque in particle 1
                torque2 -= derivativeAngluarPotential;
            }
        }
        return {force + force1, torque1, -1.0*force + force2, torque2};
    }


    /* Checks is MSM must be deactivated in certain particles. In this case, if bounded or close to bounded
     * disable MSM in particle withs type 1 and state 0. This needs to be hardcoded here for each example. */
    void patchyProteinMarkovSwitch::enableDisableMSM(vec3<double>relPosition, particleMS &part1, particleMS &part2) {
        if (relPosition.norm() <= 1.2 && part2.type == 1 && part2.state == 0) {
            part2.activeMSM = false;
        } else {
            part2.activeMSM = true;
        }
    }

    /* Given two quaternions/orientations, returns planes(unit vectors) to be aligned by torque. These
     * may vary depending on the physical arrangement of your molecules and how the patches are ordered and
     * selected. Currently set up for particle 1 with 6 binding patches and particle 2 only with one.*/
    std::tuple<vec3<double>, vec3<double>> patchyProteinMarkovSwitch::calculatePlanes(const particleMS &part1,
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

