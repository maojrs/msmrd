//
// Created by maojrs on 1/11/19.
//
#include "potentials/patchyProteinMarkovSwitch.hpp"
#include "quaternion.hpp"
#include "tools.hpp"

namespace msmrd {
    /*
     * Constructors inherited from patchyProtein class.
     */


    /* Set potentials parameters for overall potential. Consists of isotropic attractive and
     * repulsive parts plus two types of patches interactions. Needs to be redefined so parameters
     * do not depend on patchyProtein parent class. */
    void patchyProteinMarkovSwitch::setPotentialParameters() {
        // Check pacthes coordinates have unit norm
        for (auto &patch : patchesCoordinatesA) {
            if (patch.norm() != 1) {
                throw std::range_error("Patches coordinates must have norm one");
            }
        }
        for (auto &patch : patchesCoordinatesB) {
            if (patch.norm() != 1) {
                throw std::range_error("Patches coordinates must have norm one");
            }
        }

        // Set strengths of potential parts
        epsRepulsive = 1.0*strength;
        epsAttractive = -0.15*strength;
        epsPatches[0] = -0.15*strength;
        epsPatches[1] = -0.20*strength; // Special binding site

        // Set stiffness of potentials parts
        aRepulsive = 1.5;
        aAttractive = 0.75;
        aPatches[0] = 40.0;
        aPatches[1] = 40.0; // Special binding site

        // Set range parameter potentials parts
        rstarRepulsive = 0.75*sigma;
        rstarAttractive = 0.85*sigma;
        rstarPatches[0] = 0.1*sigma;
        rstarPatches[1] = 0.1*sigma; // Special binding site
    }


    /* Checks is MSM must be deactivated in certain particles. In this case, if bounded or close to bounded
     * disable MSM in particle withs type 1 and state 0. This needs to be hardcoded here for each example. */
    void patchyProteinMarkovSwitch::enableDisableMSM(particleMS &part1, particleMS &part2) {
        vec3<double> distance = part1.position - part2.position;
        if (distance.norm() <= 1.2 && part2.type == 1 && part2.state == 0) {
                part2.activeMSM = false;
        } else {
                part2.activeMSM = true;
        }
    }

    // Evaluates potential at given positions and orientations of two particles
    double patchyProteinMarkovSwitch::evaluate(particleMS &part1, particleMS &part2) {
        vec3<double> pos1 = part1.position;
        vec3<double> pos2 = part2.position;
        quaternion<double> theta1 = part1.orientation;
        quaternion<double> theta2 = part2.orientation;

        std::vector<vec3<double>> patchesCoords1;
        std::vector<vec3<double>> patchesCoords2;
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

        /* Assign patch pattern depending on particle type (note only two types of particles are supported here) */
        patchesCoords1 = assignPatches(part1.type);
        patchesCoords2 = assignPatches(part2.type);

        /* Enables/Disables MSM following potential hardcoded rules. In this case, if bounded or close to bounded
         * disable the MSM on particles type 1 and state 0. */
        enableDisableMSM(part1, part2);

        // Interaction if particles are different and if part2 in attractive state 0
        if (rvec.norm() <= 2*sigma && part2.state == 0) {
            // Loop over all patches
            for (int i = 0; i < patchesCoords1.size(); i++) {
                patchNormal1 = msmrdtools::rotateVec(patchesCoords1[i], theta1);
                patch1 = pos1 + 0.5*sigma*patchNormal1;
                for (int j = 0; j < patchesCoords2.size(); j++) {
                    patchNormal2 = msmrdtools::rotateVec(patchesCoords2[j], theta2);
                    patch2 = pos2 + 0.5*sigma*patchNormal2;
                    rpatch = patch2 - patch1; // Scale unit distance of patches by sigma
                    // Assumes the first patch from type1 has a different type of interaction,
                    if ( (i == 0 && part1.type == 0) || (j == 0 && part2.type == 0) ) {
                        patchesPotential += quadraticPotential(rpatch.norm(), sigma, epsPatches[1], aPatches[1], rstarPatches[1]);
                    }
                        // while all the other patches combinations have another type of interaction.
                    else{
                        patchesPotential += quadraticPotential(rpatch.norm(), sigma, epsPatches[0], aPatches[0], rstarPatches[0]);
                    }
                }
            }
        }

        return repulsivePotential + attractivePotential + patchesPotential;
    }

    /* Calculate and return (force1, torque1, force2, torque2), which correspond to the force and torque
     * acting on particle1 and the force and torque acting on particle2, respectively. */
    std::array<vec3<double>, 4> patchyProteinMarkovSwitch::forceTorque(particleMS &part1, particleMS &part2) {
        vec3<double> pos1 = part1.position;
        vec3<double> pos2 = part2.position;
        quaternion<double> theta1 = part1.orientation;
        quaternion<double> theta2 = part2.orientation;

        vec3<double> force;
        vec3<double> force1 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> force2 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> torque1 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> torque2 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> rvec = pos2 - pos1;
        std::vector<vec3<double>> patchesCoords1;
        std::vector<vec3<double>> patchesCoords2;
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
        // Attractive force only on conformation state 0 from particle 2; otherwise 0
        attractiveForceNorm = derivativeQuadraticPotential(rvec.norm(), sigma, epsAttractive, aAttractive, rstarAttractive);
        force = (repulsiveForceNorm + attractiveForceNorm)*rvec/rvec.norm();

        /* Enables/Disables MSM following potential hardcoded rules. In this case, if bounded or close to bounded
         * disable the MSM on particles type 1 and state 0. */
        enableDisableMSM(part1, part2);

        /* Assign patch pattern depending on particle type (note only two types of particles are supported here) */
        patchesCoords1 = assignPatches(part1.type);
        patchesCoords2 = assignPatches(part2.type);

        // Calculate forces and torque due to patches interaction if part2 in attractive state 0
        if (rvec.norm() <= 2*sigma && part2.state == 0) {
            // Loop over all patches of particle 1
            for (int i = 0; i < patchesCoords1.size(); i++) {
                patchNormal1 = msmrdtools::rotateVec(patchesCoords1[i], theta1);
                patchNormal1 = patchNormal1/patchNormal1.norm();
                patch1 = pos1 + 0.5*sigma*patchNormal1;
                // Loop over all patches of particle 2
                for (int j = 0; j < patchesCoords2.size(); j++) {
                    patchNormal2 = msmrdtools::rotateVec(patchesCoords2[j], theta2);
                    patchNormal2 = patchNormal2/patchNormal2.norm();
                    patch2 = pos2 + 0.5*sigma*patchNormal2;
                    rpatch = patch2 - patch1;
                    /* Calculate force vector between patches , correct sign of force given by rpatch/rpatch.norm().
                     * It also assumes the first patch from type = 0 has a different type of interaction. */
                    if ( (i == 0 && part1.type == 0) || (j == 0 && part2.type == 0) ) {
                        patchesForceNorm = derivativeQuadraticPotential(rpatch.norm(), sigma, epsPatches[1], aPatches[1], rstarPatches[1]);
                    }
                    else {
                        patchesForceNorm = derivativeQuadraticPotential(rpatch.norm(), sigma, epsPatches[0], aPatches[0], rstarPatches[0]);
                    }
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

    // Additional function to call forcetorque function from pybind
    std::vector<std::vector<double>>
    patchyProteinMarkovSwitch::forceTorquePyBind(particleMS &part1, particleMS &part2) {
        std::array<vec3<double>, 4> forceTorquex = forceTorque(part1, part2);
        return msmrdtools::array2Dtovec2D(forceTorquex);
    }


}