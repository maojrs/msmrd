//
// Created by maojrs on 4/29/21.
//
#include "potentials/patchyParticleSTMV.hpp"


namespace msmrd {

    patchyParticleSTMV::patchyParticleSTMV(double strength, double angularStrength)
            :  strength(strength), angularStrength(angularStrength){

        setPotentialParameters();
    };

    /* Set potentials parameters for overall potential. Consists of isotropic repulsive part
     * on m1,m2 and m3 plus three types of patches interactions*/
    void patchyParticleSTMV::setPotentialParameters() {
        // Set diameter of structure particles m1,m2 and m3:
        sigma[0] = 2 * 0.856; //m1
        sigma[1] = 2 * 1.358; //m2
        sigma[2] = 2 * 0.818; //m3

        // Set radii of all patches:
        sigmaPatches = 2 * 0.4; //i1 to i5

        // Set strengths of potential parts
        epsRepulsive = 1.0*strength;
        epsPatches[0] = -0.15*strength;
        epsPatches[1] = -0.15*strength;
        epsPatches[2] = -0.15*strength;

        // Set stiffness of potentials parts
        aRepulsive = 1.5;
        aPatches[0] = 40.0;
        aPatches[1] = 40.0;
        aPatches[2] = 40.0; // Special binding site

        // Set range parameter potentials parts
        rstarRepulsive = 0.75*sigmaPatches;
        rstarPatches[0] = 0.1*sigmaPatches;
        rstarPatches[1] = 0.1*sigmaPatches;
        rstarPatches[2] = 0.1*sigmaPatches;

        // Set coordinates of structure particles m1 to m3:
        structureCoordinates.resize(3);
        structureCoordinates[0] = vec3<double>(1.815254758, 0, 0);
        structureCoordinates[1] = vec3<double>(0, 0, 0);
        structureCoordinates[2] = vec3<double>(-0.586991122, 1.460740195, 0);

        // set coordinates of patches i1 to i5:
        patchesCoordinates.resize(5);
        patchesCoordinates[0] = vec3<double>(1.694500527, 0.674187053, -0.694570924);
        patchesCoordinates[1] = vec3<double>(2.198599972, -0.207424062, -0.826036629);
        patchesCoordinates[2] = vec3<double>(-0.487943602, -0.924886815, -0.163915576);
        patchesCoordinates[3] = vec3<double>(-1.348049652, 0.841422688, -0.057270532);
        patchesCoordinates[4] = vec3<double>(-1.003689633, 1.739415955, -0.562472092);
    }

    double patchyParticleSTMV::evaluate(particle &part1, particle &part2) {
        double repulsivePotential = 0.0;
        double patchesPotential = 0.0;

        std::array<vec3<double>, 2> relPos = relativePositionComplete(part1.position, part2.position);
        vec3<double> pos1virtual = relPos[0]; // virtual pos1 if periodic boundary; otherwise pos1.
        vec3<double> rvec = relPos[1]; //part2.position - part1.position;

        // Loop over all structure particles m1,m2 and m3 to calculate repulsive potential
        for (int i = 0; i < structureCoordinates.size(); i++) {
            auto mVecA = msmrdtools::rotateVec(structureCoordinates[i], part1.orientation);
            auto mA = pos1virtual + mVecA;
            for (int j = 0; j < structureCoordinates.size(); j++) {
                auto mVecB = msmrdtools::rotateVec(structureCoordinates[j], part2.orientation);
                auto mB = part2.position + mVecB;
                auto relM = mB - mA;
                auto sigmaTotal = 0.5 * (sigma[i]+sigma[j]);
                repulsivePotential += quadraticPotential(relM.norm(), sigmaTotal,
                        epsRepulsive, aRepulsive, rstarRepulsive);
            }
        }

        // Loop over all patches to calculate patches potential
        for (int i = 0; i < patchesCoordinates.size(); i++) {
            auto patchVec1 = msmrdtools::rotateVec(patchesCoordinates[i], part1.orientation);
            auto patch1 = pos1virtual + patchVec1;
            for (int j = 0; j < patchesCoordinates.size(); j++) {
                auto index = i+j;
                if (index== 1 or index == 5 or index == 8) {
                    auto patchVec2 = msmrdtools::rotateVec(patchesCoordinates[j], part2.orientation);
                    auto patch2 = part2.position + patchVec2;
                    auto relpatch = patch2 - patch1; // Scale unit distance of patches by sigma
                    if (index == 1){ // interaction type 1 (i1:i2)
                        patchesPotential += patchPotentialScaling * quadraticPotential(relpatch.norm(),
                                sigmaPatches,epsPatches[0], aPatches[0],rstarPatches[0]);
                    }
                    else if (index == 5){ //interaction type 2 (i3:i4)
                        patchesPotential += patchPotentialScaling * quadraticPotential(relpatch.norm(),
                                sigmaPatches,epsPatches[1], aPatches[1],rstarPatches[1]);
                    }
                    else { // interaction type 3 (i5:i5)
                        patchesPotential += patchPotentialScaling * quadraticPotential(relpatch.norm(),
                                sigmaPatches,epsPatches[2], aPatches[2],rstarPatches[2]);
                    }
                }
            }
        }
        return repulsivePotential + patchesPotential;
    };

    std::array<vec3<double>, 4> patchyParticleSTMV::forceTorque(particle &part1, particle &part2) {
        // Calculate relative position
        std::array<vec3<double>, 2> relPos = relativePositionComplete(part1.position, part2.position);
        vec3<double> pos1virtual = relPos[0]; // virtual part1.position if periodic boundary; otherwise part1.position.
        vec3<double> rvec = relPos[1]; //part2.position - part1.position;

        auto forcTorqStructure = forceTorqueStructure(part1, part2, pos1virtual);
        auto force1 = forcTorqStructure[0];
        auto torque1 = forcTorqStructure[1];
        auto force2 = forcTorqStructure[2];
        auto torque2 = forcTorqStructure[3];

        // Calculate forces and torque due to patches interaction if particles are close enough
        auto forcTorqPatches = forceTorquePatches(part1, part2, pos1virtual);
        auto force1patches = forcTorqPatches[0];
        auto torque1patches = forcTorqPatches[1];
        auto force2patches = forcTorqPatches[2];
        auto torque2patches = forcTorqPatches[3];
        return {force1 + force1patches, torque1 + torque1patches,
                force2 + force2patches, torque2 + torque2patches};
    };


    /* Calculates forces and torques due to structural interactions, first two vectors returned are the force
     * and torque applied to particle 1 and the second two vectors are the force and torque applied to particle 2.
     * This function is called by main forceTorque function. */
    std::array<vec3<double>, 4> patchyParticleSTMV::forceTorqueStructure(
            particle &part1, particle &part2, const vec3<double> pos1virtual) {

        vec3<double> repulsiveForce;
        double repulsiveForceNorm = 0.0;

        vec3<double> force1 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> force2 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> torque1 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> torque2 = vec3<double> (0.0, 0.0, 0.0);

        // Loop over all structure particles m1,m2 and m3 to calculate repulsive potential
        for (int i = 0; i < structureCoordinates.size(); i++) {
            auto mVecA = msmrdtools::rotateVec(structureCoordinates[i], part1.orientation);
            auto mA = pos1virtual + mVecA;
            // Loop over all patches of particle 2
            for (int j = 0; j < structureCoordinates.size(); j++) {
                auto mVecB = msmrdtools::rotateVec(structureCoordinates[j], part2.orientation);
                auto mB = part2.position + mVecB;
                auto relM = mB - mA;
                auto sigmaTotal = 0.5 * (sigma[i]+sigma[j]);

                // Calculate force vector between patches , correct sign of force given by relpatch/relpatch.norm().
                repulsiveForceNorm = derivativeQuadraticPotential(relM.norm(),
                        sigmaTotal, epsRepulsive, aRepulsive, rstarRepulsive);

                if (relM.norm() == 0) {
                    repulsiveForce = vec3<double>(0, 0, 0);
                } else {
                    repulsiveForce = repulsiveForceNorm * relM / relM.norm();
                }
                // Calculate force and torque acting on particle 1 and add values to previous forces and torques
                force1 += patchPotentialScaling * repulsiveForce;
                torque1 += patchPotentialScaling * mVecA.cross(repulsiveForce);

                // Calculate force and torque acting on particle 2 and add values to previous forces and torques
                force2 += -1.0 * patchPotentialScaling * repulsiveForce;
                torque2 += patchPotentialScaling * mVecB.cross(-1.0 * repulsiveForce);
            }
        }
        return {force1, torque1, force2, torque2};
    };


    /* Calculates forces and torques due to pacthes interactions, first two vectors returned are the force
     * and torque applied to particle 1 and the second two vectors are the force and torque applied to particle 2.
     * This function is called by main forceTorque function. */
    std::array<vec3<double>, 4> patchyParticleSTMV::forceTorquePatches(
            particle &part1, particle &part2, const vec3<double> pos1virtual) {

        vec3<double> patchForce;
        double patchesForceNorm = 0.0;

        vec3<double> force1 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> force2 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> torque1 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> torque2 = vec3<double> (0.0, 0.0, 0.0);

        // Loop over all patches of particle 1
        for (int i = 0; i < patchesCoordinates.size(); i++) {
            auto patchVec1 = msmrdtools::rotateVec(patchesCoordinates[i], part1.orientation);
            auto patch1 = pos1virtual + patchVec1;
            // Loop over all patches of particle 2
            for (int j = 0; j < patchesCoordinates.size(); j++) {
                auto index = i + j;
                if (index == 1 or index == 5 or index == 8) {
                    auto patchVec2 = msmrdtools::rotateVec(patchesCoordinates[j], part2.orientation);
                    auto patch2 = part2.position + patchVec2;
                    auto relpatch = patch2 - patch1;
                    if (index == 1) { // interaction type 1 (i1:i2)
                        // Calculate force vector between patches , correct sign of force given by relpatch/relpatch.norm().
                        patchesForceNorm = derivativeQuadraticPotential(relpatch.norm(),
                                sigmaPatches, epsPatches[0], aPatches[0], rstarPatches[0]);
                    }
                    else if (index == 5) { // interaction type 2 (i3:i4)
                        patchesForceNorm = derivativeQuadraticPotential(relpatch.norm(),
                                sigmaPatches, epsPatches[1], aPatches[1], rstarPatches[1]);
                    }
                    else{ // interaction type 2 (i5:i5)
                        patchesForceNorm = derivativeQuadraticPotential(relpatch.norm(),
                                sigmaPatches, epsPatches[2], aPatches[2], rstarPatches[2]);
                    }
                    if (relpatch.norm() == 0) {
                        patchForce = vec3<double>(0, 0, 0);
                    } else {
                        patchForce = patchesForceNorm * relpatch / relpatch.norm();
                    }
                    // Calculate force and torque acting on particle 1 and add values to previous forces and torques
                    force1 += patchPotentialScaling * patchForce;
                    torque1 += patchPotentialScaling * patchVec1.cross(patchForce);

                    // Calculate force and torque acting on particle 2 and add values to previous forces and torques
                    force2 += -1.0 * patchPotentialScaling * patchForce;
                    torque2 += patchPotentialScaling * patchVec2.cross(-1.0 * patchForce);
                }
            }
        }
        return {force1, torque1, force2, torque2};
    };

    /* Custom quadratic potential (same functions as in patchy particle, inheritance was not an option :( ):
     * @param r is distance between particles or  patches,
     * @param eps strength of the potential
     * @param rstar distance to switch to second potential function
     * @param a the stiffness, together with rstar determines the range of the potential */
    double patchyParticleSTMV::quadraticPotential(double r, double sig, double eps, double a, double rstar) {
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

    // Derivative of patchyParticleSTMV::quadraticPotential
    double patchyParticleSTMV::derivativeQuadraticPotential(double r, double sig, double eps, double a, double rstar) {
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

}
