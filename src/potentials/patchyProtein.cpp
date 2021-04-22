//
// Created by maojrs on 9/26/18.
//

#include <utility>
#include "potentials/patchyProtein.hpp"
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
    patchyProtein::patchyProtein(double sigma, double strength,
                                 std::vector<vec3<double>> patchesCoordinatesA,
                                 std::vector<vec3<double>> patchesCoordinatesB)
            :  sigma(sigma), strength(strength),
               patchesCoordinatesA(std::move(patchesCoordinatesA)),
               patchesCoordinatesB(std::move(patchesCoordinatesB)) {

       setPotentialParameters();
    };

    patchyProtein::patchyProtein(double sigma, double strength,
                                 std::vector<std::vector<double>> patchesCoordsA,
                                 std::vector<std::vector<double>> patchesCoordsB)
            :  sigma(sigma), strength(strength) {
        patchesCoordinatesA.resize(patchesCoordsA.size());
        patchesCoordinatesB.resize(patchesCoordsB.size());
        for (int i=0; i < patchesCoordsA.size(); i++ ) {
            patchesCoordinatesA[i] = vec3<double> (patchesCoordsA[i]);
        }
        for (int i=0; i < patchesCoordsB.size(); i++ ) {
            patchesCoordinatesB[i] = vec3<double> (patchesCoordsB[i]);
        }
        setPotentialParameters();
    }


    /* Set potentials parameters for overall potential. Consists of isotropic attractive and
     * repulsive parts plus two types of patches interactions*/
    void patchyProtein::setPotentialParameters() {
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

        // Check there are patches, if not, turn patchy interaction off
        if (patchesCoordinatesA.size() == 0 or patchesCoordinatesB.size() == 0) {
            patchesActive = false;
        }
    }


    // Evaluates potential at given positions and orientations of two particles
    double patchyProtein::evaluate(particle &part1, particle &part2) {

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

        // Evaluate patches potential if particles are close enough
        double patchesPotential = 0.0;
        if (rvec.norm() <= 2*sigma and patchesActive) {
            // Evaluate patches potential, using auxiliary function
            patchesPotential = evaluatePatchesPotential(part1, part2, pos1virtual, patchesCoords1, patchesCoords2);
        }

        return repulsivePotential + attractivePotential + patchesPotential;

    }


    /* Auxiliary function that calculates patches interaction contribution to potential. Called by main
     * evaluate function. */
    double patchyProtein::evaluatePatchesPotential(particle &part1, particle &part2,
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

        // Loop over all patches
        for (int i = 0; i < patchesCoords1.size(); i++) {
            patchNormal1 = msmrdtools::rotateVec(patchesCoords1[i], part1.orientation);
            patch1 = pos1virtual + 0.5 * sigma * patchNormal1;
            for (int j = 0; j < patchesCoords2.size(); j++) {
                patchNormal2 = msmrdtools::rotateVec(patchesCoords2[j], part2.orientation);
                patch2 = part2.position + 0.5 * sigma * patchNormal2;
                rpatch = patch2 - patch1; // Scale unit distance of patches by sigma
                // Assumes the first patch from type "0" has a different type of interaction,
                if ((i == 0 && part1.type == 0) || (j == 0 && part2.type == 0)) {
                    patchesPotential += quadraticPotential(rpatch.norm(), sigma, epsPatches[1],
                                                           aPatches[1], rstarPatches[1]);
                }
                // while all the other patches combinations have another type of interaction.
                else {
                    patchesPotential += quadraticPotential(rpatch.norm(), sigma, epsPatches[0],
                                                           aPatches[0], rstarPatches[0]);
                }
            }
        }

        return patchPotentialScaling * patchesPotential;
    };


    /* Calculate and return (force1, torque1, force2, torque2), which correspond to the force and torque
     * acting on particle1 and the force and torque acting on particle2, respectively. */
    std::array<vec3<double>, 4> patchyProtein::forceTorque(particle &part1, particle &part2) {

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

        // Calculate forces and torque due to patches interaction if particles are close enough
        if (rvec.norm() <= 2*sigma and patchesActive) {
            auto forcTorqPatches = forceTorquePatches(part1, part2, pos1virtual, patchesCoords1, patchesCoords2);
            auto force1 = forcTorqPatches[0];
            auto torque1 = forcTorqPatches[1];
            auto force2 = forcTorqPatches[2];
            auto torque2 = forcTorqPatches[3];
            return {force + force1, torque1, -1.0 * force + force2, torque2};
        } else {
            // Return only isotropic force if particles are far
            auto zeroTorque = vec3<double> (0.0, 0.0, 0.0);
            return {force, zeroTorque, -1.0 * force, zeroTorque};
        }
    }

    /* Auxiliary function that calculates patches interaction forces. Called by main forceTorque function. */
    std::array<vec3<double>, 4> patchyProtein::forceTorquePatches(particle &part1, particle &part2,
                                                                  vec3<double> &pos1virtual,
                                                                  std::vector<vec3<double>> &patchesCoords1,
                                                                  std::vector<vec3<double>> &patchesCoords2){
        // Initialize forceas and torques due to patches
        vec3<double> force1 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> force2 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> torque1 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> torque2 = vec3<double> (0.0, 0.0, 0.0);

        // auxiliary variables to calculate force and torque
        double patchesForceNorm = 0.0;
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
                 * It also assumes the first patch from type = 0 has a different type of interaction. */
                if ((i == 0 && part1.type == 0) || (j == 0 && part2.type == 0)) {
                    patchesForceNorm = derivativeQuadraticPotential(rpatch.norm(), sigma, epsPatches[1],
                                                                    aPatches[1], rstarPatches[1]);
                } else {
                    patchesForceNorm = derivativeQuadraticPotential(rpatch.norm(), sigma, epsPatches[0],
                                                                    aPatches[0], rstarPatches[0]);
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
    };


    /* Custom quadratic potential (same functions as in patchy particle, inheritance was not an option :( ):
     * @param r is distance between particles or  patches,
     * @param eps strength of the potential
     * @param rstar distance to switch to second potential function
     * @param a the stiffness, together with rstar determines the range of the potential */
    double patchyProtein::quadraticPotential(double r, double sig, double eps, double a, double rstar) {
        // Parameters to force continuity and continuous derivative
        double rcritical = std::pow(sig, 2)/(a*rstar);
        double b = (1.0 - a*std::pow(rstar/sig, 2))/std::pow(rcritical/sig - rstar/sig, 2);
        // Distance cannot be negative, so it is always the absolute value
        double rlocal = std::abs(r);
        // Calculate different regions of potential
        if (rlocal < rstar) {
            return eps * (1.0 - a * std::pow((rlocal / sig),2));
        } else if (rlocal < rcritical) {
            return eps * b * std::pow((rcritical / sig - rlocal / sig), 2);
        } else {
            return 0.0;
        }
    }

    // Derivative of patchyProtein::quadraticPotential
    double patchyProtein::derivativeQuadraticPotential(double r, double sig, double eps, double a, double rstar) {
        // Parameters to force continuity and continuous derivative
        double rcritical = std::pow(sig, 2)/(a*rstar);
        double b = (1.0 - a*std::pow(rstar/sig, 2))/std::pow(rcritical/sig - rstar/sig, 2);
        // Distance cannot be negative, so it is always the absolute value
        double rlocal = std::abs(r);
        // Calculate derivative in different regions of potential
        if (rlocal < rstar) {
            return - 2.0*eps*a * rlocal/std::pow(sig,2);
        } else if (rlocal < rcritical) {
            return -2.0*eps * b * (rcritical - rlocal)/std::pow(sig,2);
        } else {
            return 0.0;
        }
    }

    // Assign pacthes coordinates in terms of particle type
    std::vector<vec3<double>> patchyProtein::assignPatches(int type) {
        if (type == 0) {
            return patchesCoordinatesA;
        }
        else if (type == 1) {
            return patchesCoordinatesB;
        }
        else {
            throw std::runtime_error("This potential only supports two different types of particles.");
        }
    }


}
