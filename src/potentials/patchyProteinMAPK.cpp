//
// Created by maojrs on 3/18/21.
//


#include <utility>
#include "potentials/patchyProteinMAPK.hpp"
#include "quaternion.hpp"
#include "tools.hpp"

namespace msmrd {

    /*
    * Constructors analogous to the patchyParticle potential.
    * @param sigma diameter of sphere at which patches are placed. Corresponds to the radius of the particle.
    * @param strength overall strength of the potential
    * @param patchesCoordintatesA list of patches coordinates in a unit sphere corresponding to type1 particles.
    * @param patchesCoordintatesB list of patches coordinates in a unit sphere corresponding to type2 particles

    */
    patchyProteinMAPK::patchyProteinMAPK(double sigma, double strength,
                                 std::vector<vec3<double>> patchesCoordinatesA,
                                 std::vector<vec3<double>> patchesCoordinatesB)
            :  patchyProtein(sigma, strength, patchesCoordinatesA, patchesCoordinatesB ) {

        setPotentialParameters();
    };

    patchyProteinMAPK::patchyProteinMAPK(double sigma, double strength,
                                 std::vector<std::vector<double>> patchesCoordsA,
                                 std::vector<std::vector<double>> patchesCoordsB)
            :  patchyProtein(sigma, strength, patchesCoordsA, patchesCoordsB ) {
        setPotentialParameters();
    }

    /* Set potentials parameters for overall potential. Consists of isotropic attractive and
     * repulsive parts plus two types of patches interactions*/
    void patchyProteinMAPK::setPotentialParameters() {
        // Check pacthes coordinates have unit norm
        for (auto &patch : patchesCoordinatesA) {
            if (patch.norm() != 1) {
                throw std::invalid_argument("Patches coordinates must have norm one");
            }
        }
        for (auto &patch : patchesCoordinatesB) {
            if (patch.norm() != 1) {
                throw std::invalid_argument("Patches coordinates must have norm one");
            }
        }

        // Set strengths of potential parts
        epsRepulsive = 1.0*strength;
        epsAttractive = 0.0; //-0.15*strength;
        epsPatches[0] = -0.15*strength; // When interacting with Kinase
        epsPatches[1] = -0.15*strength; // When interacting with Phosphotase

        // Set stiffness of potentials parts
        aRepulsive = 1.5;
        aAttractive = 0.75;
        aPatches[0] = 40.0; // When interacting with Kinase
        aPatches[1] = 40.0; // When interacting with Phosphotase

        // Set range parameter potentials parts
        rstarRepulsive = 0.75*sigma;
        rstarAttractive = 0.85*sigma;
        rstarPatches[0] = 0.1*sigma; // When interacting with Kinase
        rstarPatches[1] = 0.1*sigma; // When interacting with Phosphotase

        // Check there are patches, if not, turn patchy interaction off
        if (patchesCoordinatesA.size() == 0 or patchesCoordinatesB.size() == 0) {
            patchesActive = false;
        }
    }

    /* Auxiliary function that calculates patches interaction contribution to potential. Called by main
     * evaluate function. */
    double patchyProteinMAPK::evaluatePatchesPotential(particle &part1, particle &part2,
                                                   vec3<double> &pos1virtual,
                                                   std::vector<vec3<double>> &patchesCoords1,
                                                   std::vector<vec3<double>> &patchesCoords2){
        // Declare variables used in loop
        double patchesPotential = 0.0;
        vec3<double> patchNormal1;
        vec3<double> patchNormal2;
        vec3<double> patch1;
        vec3<double> patch2;
        vec3<double> rpatch;
        bool interact = false;

        // Loop over all patches
        for (int i = 0; i < patchesCoords1.size(); i++) {
            patchNormal1 = msmrdtools::rotateVec(patchesCoords1[i], part1.orientation);
            patch1 = pos1virtual + 0.5 * sigma * patchNormal1;
            for (int j = 0; j < patchesCoords2.size(); j++) {
                patchNormal2 = msmrdtools::rotateVec(patchesCoords2[j], part2.orientation);
                patch2 = part2.position + 0.5 * sigma * patchNormal2;
                rpatch = patch2 - patch1; // Scale unit distance of patches by sigma
                /* Only allow for interactions between MAPK and kinase or MAPK and phosphosate
                 * Different type of interaction with kinase, state 0 corresponds to active state. */
                // part1 is MAPK and part2 is kinase:
                if (part1.type == 0 and part2.type == 1 and part2.state == 0) {
                    if ((i == 0 and (part1.state == 0 or part1.state == 2)) or
                        (i == 1 and (part1.state == 0 or part1.state == 1))) {
                        patchesPotential += quadraticPotential(rpatch.norm(), sigma, epsPatches[0],
                                                               aPatches[0], rstarPatches[0]);
                    }
                }
                // part2 is MAPK and part1 is kinase:
                else if (part2.type == 0 and part1.type == 1 and part1.state == 0) {
                    if ((j == 0 and (part2.state == 0 or part2.state == 2)) or
                        (j == 1 and (part2.state == 0 or part2.state == 1))) {
                        patchesPotential += quadraticPotential(rpatch.norm(), sigma, epsPatches[0],
                                                               aPatches[0], rstarPatches[0]);
                    }
                }
                // part1 is MAPK and part2 is phosphotase:
                else if (part1.type == 0  and part2.type == 2 and part2.state == 0) {
                    if ((i == 0 and (part1.state == 1 or part1.state == 4)) or
                        (i == 1 and (part1.state == 2 or part1.state == 4))) {
                        patchesPotential += quadraticPotential(rpatch.norm(), sigma, epsPatches[1],
                                                               aPatches[1], rstarPatches[1]);
                    }
                }
                // part2 is MAPK and part1 is phosphotase:
                else if (part2.type == 0  and part1.type == 2 and part1.state == 0) {
                    if ((j == 0 and (part2.state == 1 or part2.state == 4)) or
                        (j == 1 and (part2.state == 2 or part2.state == 4))) {
                        patchesPotential += quadraticPotential(rpatch.norm(), sigma, epsPatches[1],
                                                               aPatches[1], rstarPatches[1]);
                    }
                }
            }
        }

        return patchPotentialScaling * patchesPotential;
    }


    /* Auxiliary function that calculates patches interaction forces. Called by main forceTorque function. */
    std::array<vec3<double>, 4> patchyProteinMAPK::forceTorquePatches(particle &part1, particle &part2,
                                                                  vec3<double> &pos1virtual,
                                                                  std::vector<vec3<double>> &patchesCoords1,
                                                                  std::vector<vec3<double>> &patchesCoords2){
        // Initialize forceas and torques due to patches
        vec3<double> force1 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> force2 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> torque1 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> torque2 = vec3<double> (0.0, 0.0, 0.0);

        // auxiliary variables to calculate force and torque
        double patchesForceNorm;
        vec3<double> patchForce;
        vec3<double> patch1;
        vec3<double> patch2;
        vec3<double> patchNormal1;
        vec3<double> patchNormal2;
        vec3<double> rpatch;

        // Loop over all patches of particle 1
        for (int i = 0; i < patchesCoords1.size(); i++) {
            patchNormal1 = msmrdtools::rotateVec(patchesCoords1[i], part1.orientation);
            patchNormal1 = patchNormal1 / patchNormal1.norm();
            patch1 = pos1virtual + 0.5 * sigma * patchNormal1;
            // Loop over all patches of particle 2
            for (int j = 0; j < patchesCoords2.size(); j++) {
                patchNormal2 = msmrdtools::rotateVec(patchesCoords2[j], part2.orientation);
                patchNormal2 = patchNormal2 / patchNormal2.norm();
                patch2 = part2.position + 0.5 * sigma * patchNormal2;
                // Calculate distance between the two patches
                rpatch = patch2 - patch1;
                /* Calculate force vector between patches , correct sign of force given by rpatch/rpatch.norm().
                 * Only allow for interactions between MAPK and kinase or MAPK and phosphosate .It also assumes
                 * different type of interaction with kinase/phosphotase, and state 0 corresponds to active state. */
                // part1 is MAPK and part2 is kinase:
                if (part1.type == 0 and part2.type == 1 and part2.state == 0) {
                    if ((i == 0 and (part1.state == 0 or part1.state == 2)) or
                        (i == 1 and (part1.state == 0 or part1.state == 1))) {
                        patchesForceNorm = derivativeQuadraticPotential(rpatch.norm(), sigma, epsPatches[0],
                                                                        aPatches[0], rstarPatches[0]);
                    }
                }
                // part2 is MAPK and part1 is kinase:
                else if (part2.type == 0 and part1.type == 1 and part1.state == 0) {
                    if ((j == 0 and (part2.state == 0 or part2.state == 2)) or
                        (j == 1 and (part2.state == 0 or part2.state == 1))) {
                        patchesForceNorm = derivativeQuadraticPotential(rpatch.norm(), sigma, epsPatches[0],
                                                                        aPatches[0], rstarPatches[0]);
                    }
                }
                // part1 is MAPK and part2 is phosphotase:
                else if (part1.type == 0  and part2.type == 2 and part2.state == 0) {
                    if ((i == 0 and (part1.state == 1 or part1.state == 3)) or
                        (i == 1 and (part1.state == 2 or part1.state == 3))) {
                        patchesForceNorm = derivativeQuadraticPotential(rpatch.norm(), sigma, epsPatches[1],
                                                                        aPatches[1], rstarPatches[1]);
                    }
                }
                // part2 is MAPK and part1 is phosphotase:
                else if (part2.type == 0  and part1.type == 2 and part1.state == 0) {
                    if ((j == 0 and (part2.state == 1 or part2.state == 3)) or
                        (j == 1 and (part2.state == 2 or part2.state == 3))) {
                        patchesForceNorm = derivativeQuadraticPotential(rpatch.norm(), sigma, epsPatches[1],
                                                                        aPatches[1], rstarPatches[1]);
                    }
                }
                // Determine force vector avoiding division by zero
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

        return {force1, torque1, force2, torque2};
    }

    // Assign pacthes coordinates in terms of particle type
    std::vector<vec3<double>> patchyProteinMAPK::assignPatches(int type) {
        if (type == 0) {
            return patchesCoordinatesA;
        } else if (type == 1 or type == 2) {
            return patchesCoordinatesB;
        } else {
            throw std::runtime_error("This potential only supports three different types of particles.");
        }
    }


}

