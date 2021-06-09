//
// Created by maojrs on 4/29/21.
//
#include "potentials/patchyParticleSTMV.hpp"


namespace msmrd {

    patchyParticleSTMV::patchyParticleSTMV(double strength, double angularStrength)
            :  strength(strength), angularStrength(angularStrength){

        setPotentialParameters();
        setBoundStates();
    };

    /* Set potentials parameters for overall potential. Consists of isotropic repulsive part
     * on m1,m2 and m3 plus three types of patches interactions*/
    void patchyParticleSTMV::setPotentialParameters() {
        // Set diameter of structure particles m1,m2 and m3:
        sigma[0] = 2 * 0.856; //m1 2.2 to move binding site slightly away from patch by 0.4 (0.2 on each side)
        sigma[1] = 2 * 1.358; //m2
        sigma[2] = 2 * 0.818; //m3 2.2 to move binding site slightly away from patch by 0.4 (0.2 on each side)

        // Set diameter of all patches:
        sigmaPatches = 2 * 0.4; //i1 to i5

        // Set strengths of potential parts
        epsRepulsive = 1.0*strength;
        epsPatchesList[0] = -0.50*strength;
        epsPatchesList[1] = -1.50*strength;
        epsPatchesList[2] = -1.50*strength;

        // Set stiffness of potentials parts
        aRepulsive = 1.5;
        aPatchesList[0] = 10.0; //40.0;
        aPatchesList[1] = 10.0; //40.0;
        aPatchesList[2] = 10.0; //40.0;

        /* Set range parameter potentials parts  */
        rstarRepulsive = 0.75; // will be multiplied by corresponding sigma in implementation
        rstarPatchesList[0] = 0.15 * sigmaPatches;
        rstarPatchesList[1] = 0.15 * sigmaPatches;
        rstarPatchesList[2] = 0.15 * sigmaPatches;

        // Set coordinates of structure particles m1 to m3:
        structureCoordinates.resize(3);
        structureCoordinates[0] = vec3<double>(1.815254758, 0, 0);
        structureCoordinates[1] = vec3<double>(0, 0, 0);
        structureCoordinates[2] = vec3<double>(-0.586991122, 1.460740195, 0);

        // set coordinates of patches i1 to i5:
        patchesCoordinates.resize(5);
        patchesCoordinates[0] = 1.0*vec3<double>(1.694500527, 0.674187053, -0.694570924);
        patchesCoordinates[1] = 1.0*vec3<double>(2.198599972, -0.207424062, -0.826036629);
        patchesCoordinates[2] = 1.0*vec3<double>(-0.487943602, -0.924886815, -0.163915576);
        patchesCoordinates[3] = 1.0*vec3<double>(-1.348049652, 0.841422688, -0.057270532);
        patchesCoordinates[4] = 1.0*vec3<double>(-1.003689633, 1.739415955, -0.562472092);
    }

    void patchyParticleSTMV::setBoundStates(){
        /* Define relative position vectors from particle 1 at the origin (nm). */
        std::array<vec3<double>, 5> relPos;
        relPos[0] = { 2.39531987, -2.43096797, -1.21785841};
        relPos[1] = -1.0 * relPos[0];
        relPos[2] = {-2.13123298, -2.18629724,  0.03285758};
        relPos[3] = -1.0 * relPos[2];
        relPos[4] = {-2.06269526,  3.99529612, -1.3770492};
        /* Relative orientations */
        std::array<quaternion<double>, 6> quatRotations;
        quatRotations[0] = {0.809015  , 0.34173024, 0.10274882, 0.46707371};
        quatRotations[1] = quatRotations[0].conj();
        quatRotations[2] = {0.50000148,  0.00507126, -0.01793432, -0.86582398};
        quatRotations[3] = quatRotations[1].conj();
        quatRotations[4] = {1.27370881e-05, -8.61936819e-02, -3.64128881e-01, -9.27351501e-01};
        /* Fill bound states with corresponding combinations of relative position vectors and quaternion orientations.
         * Note we define relativeOrientation as q2 * q1.conj(), so it matches the relative orientation from
         * particle 1 as used in trajectories/discrete/discreteTrajectory.hpp */
        boundStates[0] = std::make_tuple(relPos[0], quatRotations[0]); // part1:i2, part2:i1
        boundStates[1] = std::make_tuple(relPos[1], quatRotations[1]); // part1:i1, part2:i2
        boundStates[2] = std::make_tuple(relPos[2], quatRotations[2]); // part1:i3, part2:i4
        boundStates[3] = std::make_tuple(relPos[3], quatRotations[3]); // part1:i4, part2:i3
        boundStates[4] = std::make_tuple(relPos[4], quatRotations[4]); // part1:i5, part2:i5
        /* Define two reference planes for each bound state. This reference planes are the ones
         * to be aligned by the potential. Plane one corresponds to the vector part2:ik-part1:il
         * (dihedral angle) and the plane2 is (part2:ik-part1:il).cross(part2:ik)*/
        std::array<vec3<double>, 5> plane1;
        std::array<vec3<double>, 5> plane2;
        std::vector<vec3<double>> part2PatchesCoordinates(5);
        // Bound state 0
        auto vec_i1 = relPos[0] + msmrdtools::rotateVec(patchesCoordinates[0], quatRotations[0]); // part2:i1
        plane1[0] = vec_i1 - patchesCoordinates[1];
        plane2[0] = plane1[0].cross(vec_i1);
        // Bound state 1
        auto vec_i2 = relPos[1] + msmrdtools::rotateVec(patchesCoordinates[1], quatRotations[1]); // part2:i2
        plane1[1] = vec_i2 - patchesCoordinates[0];
        plane2[1] = plane1[1].cross(vec_i2);
        // Bound state 2
        auto vec_i4 = relPos[2] + msmrdtools::rotateVec(patchesCoordinates[3], quatRotations[2]); // part2:i4
        plane1[2] = vec_i4 - patchesCoordinates[2];
        plane2[2] = plane1[2].cross(vec_i4);
        // Bound state 3
        auto vec_i3 = relPos[3] + msmrdtools::rotateVec(patchesCoordinates[2], quatRotations[3]); // part2:i3
        plane1[3] = vec_i3 - patchesCoordinates[3];
        plane2[3] = plane1[3].cross(vec_i3);
        // Bound state 4
        auto vec_i5 = relPos[4] + msmrdtools::rotateVec(patchesCoordinates[4], quatRotations[4]); // part2:i5
        plane1[4] = vec_i5 - patchesCoordinates[4];
        plane2[4] = plane1[4].cross(vec_i5);
        // Normalize and save reference planes
        for (int i=0; i<5; i++) {
            plane1[i] = plane1[i]/plane1[i].norm();
            plane2[i] = plane2[i]/plane2[i].norm();
            refPlanes[i] = std::make_tuple(plane1[i], plane2[i]);
        }
    };

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
                auto rstarTotal = sigmaTotal * rstarRepulsive;
                repulsivePotential += quadraticPotential(relM.norm(), sigmaTotal,
                        epsRepulsive, aRepulsive, rstarTotal);
            }
        }

        // If particles are close enough, calculate patches potential
        if (rvec.norm() <= minimumR) {
            patchesPotential = evaluatePatchesPotential(part1, part2, pos1virtual);
        }
        return repulsivePotential + patchesPotential;
    };

    // Just evaluate patches potential part of potential (including orientational part)
    double patchyParticleSTMV::evaluatePatchesPotential( particle &part1, particle &part2, vec3<double> &pos1virtual) {
        double patchesPotential = 0.0;
        double angularPotential = 0.0;

        // Evaluate pair potentials (no orientations)
        for (int i=0; i < patchesCoordinates.size(); i++) {
            auto patchVec1 = msmrdtools::rotateVec(patchesCoordinates[i], part1.orientation);
            auto patch1 = pos1virtual + 1.1 * patchVec1; // 1.1 factor to keep the bnding site slightly off the patch
            for (int j=1; j < patchesCoordinates.size(); j++) {
                auto patchVec2 = msmrdtools::rotateVec(patchesCoordinates[j], part2.orientation);
                auto patch2 = part2.position + 1.1 * patchVec2;
                auto relpatch = patch2 - patch1;
                if ((i == 0 and j == 1) or (i == 1 and j == 0)) { // interaction type 1 (i1:i2)
                    patchesPotential += patchPotentialScaling * quadraticPotential(relpatch.norm(),
                            sigmaPatches, epsPatchesList[0],aPatchesList[0], rstarPatchesList[0]);
                } else if ((i == 2 and j == 3) or (i == 3 and j == 2)) { //interaction type 2 (i3:i4)
                    patchesPotential += patchPotentialScaling * quadraticPotential(relpatch.norm(),
                            sigmaPatches, epsPatchesList[1],aPatchesList[1], rstarPatchesList[1]);
                } else if (i == 4 and j == 4) { // interaction type 3 (i5:i5)
                    patchesPotential += patchPotentialScaling * quadraticPotential(relpatch.norm(),
                            sigmaPatches, epsPatchesList[2],aPatchesList[2], rstarPatchesList[2]);
                }
            }
        }

        /* Evaluate orientational contributions first. If one uses only one minima, use potential
         * of -(1/16)*(cos(theta) + 1)^4 */
        vec3<double> plane1;
        vec3<double> plane2;
        vec3<double> refPlane1;
        vec3<double> refPlane2;
        int boundStateIndex;
        std::tie(plane1, plane2, boundStateIndex) = calculatePlanes(part1, part2, pos1virtual);
        if (boundStateIndex > -1) {
            std::tie(refPlane1, refPlane2) = refPlanes[boundStateIndex];
            // Potential to align plane1:
            double cosSquared = (plane1 * refPlane1 + 1) * (plane1 * refPlane1 + 1);
            angularPotential += -(1.0 / 16.0) * angularStrength * cosSquared * cosSquared;
            // Potential to align plane2:
            cosSquared = (plane2 * refPlane2 + 1) * (plane1 * refPlane2 + 1);
            angularPotential += -(1.0 / 16.0) * angularStrength * cosSquared * cosSquared;
        }
        return patchesPotential + angularPotential;
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
        if (rvec.norm() <= minimumR) {
            auto forcTorqPatches = forceTorquePatches(part1, part2, pos1virtual);
            auto force1patches = forcTorqPatches[0];
            auto torque1patches = forcTorqPatches[1];
            auto force2patches = forcTorqPatches[2];
            auto torque2patches = forcTorqPatches[3];
            return {force1 + force1patches, torque1 + torque1patches,
                    force2 + force2patches, torque2 + torque2patches};
        } else{
            return {force1, torque1, force2, torque2};
        }
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
        // Loop over all structures of particle 1
        for (int i = 0; i < structureCoordinates.size(); i++) {
            auto mVecA = msmrdtools::rotateVec(structureCoordinates[i], part1.orientation);
            auto mA = pos1virtual + mVecA;
            // Loop over all structures of particle 2
            for (int j = 0; j < structureCoordinates.size(); j++) {
                auto mVecB = msmrdtools::rotateVec(structureCoordinates[j], part2.orientation);
                auto mB = part2.position + mVecB;
                auto relM = mB - mA;
                auto sigmaTotal = 0.5 * (sigma[i]+sigma[j]);
                auto rstarTotal = sigmaTotal * rstarRepulsive;

                // Calculate force vector between patches , correct sign of force given by relpatch/relpatch.norm().
                repulsiveForceNorm = derivativeQuadraticPotential(relM.norm(),
                        sigmaTotal, epsRepulsive, aRepulsive, rstarTotal);

                if (relM.norm() == 0) {
                    repulsiveForce = vec3<double>(0, 0, 0);
                } else {
                    repulsiveForce = repulsiveForceNorm * relM / relM.norm();
                }
                // Calculate force and torque acting on particle 1 and add values to previous forces and torques
                force1 += repulsiveForce;
                torque1 += mVecA.cross(repulsiveForce);

                // Calculate force and torque acting on particle 2 and add values to previous forces and torques
                force2 += -1.0 * repulsiveForce;
                torque2 += mVecB.cross(-1.0 * repulsiveForce);
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
        double patchesForceNorm;

        vec3<double> force1 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> force2 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> torque1 = vec3<double> (0.0, 0.0, 0.0);
        vec3<double> torque2 = vec3<double> (0.0, 0.0, 0.0);

        // Loop over patches (all if no active bound only one if within bound tolerance)
        for (int i=0; i < patchesCoordinates.size(); i++) {
            auto patchVec1 = msmrdtools::rotateVec(patchesCoordinates[i], part1.orientation);
            auto patch1 = pos1virtual + 1.1 * patchVec1; // 1.1 factor to keep the bnding site slightly off the patch
            // Loop over all patches of particle 2
            for (int j=0; j < patchesCoordinates.size(); j++) {
                patchesForceNorm = 0.0;
                auto patchVec2 = msmrdtools::rotateVec(patchesCoordinates[j], part2.orientation);
                auto patch2 = part2.position + 1.1 * patchVec2; // 1.1 factor to keep the bnding site slightly off the patch
                auto relpatch = patch2 - patch1;
                if ((i == 0 and j == 1) or (i == 1 and j == 0)) { // interaction type 1 (i1:i2)
                    // Calculate force vector between patches , correct sign of force given by relpatch/relpatch.norm().
                    patchesForceNorm = derivativeQuadraticPotential(relpatch.norm(),
                            sigmaPatches, epsPatchesList[0], aPatchesList[0], rstarPatchesList[0]);
                }
                else if ((i == 2 and j == 3) or (i == 3 and j == 2)) { // interaction type 2 (i3:i4)
                    patchesForceNorm = derivativeQuadraticPotential(relpatch.norm(),
                            sigmaPatches, epsPatchesList[1], aPatchesList[1], rstarPatchesList[1]);
                }
                else if (i == 4 and j == 4){ // interaction type 2 (i5:i5)
                    patchesForceNorm = derivativeQuadraticPotential(relpatch.norm(),
                            sigmaPatches, epsPatchesList[2], aPatchesList[2], rstarPatchesList[2]);
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

        /* Next calculate orientational contribution. If one uses only one minima, use potential
         * of -(1/16)*(cos(theta) + 1)^4. It is only active if patches are closer than
         * angularPotentialTolerance */
        vec3<double> plane1;
        vec3<double> plane2;
        vec3<double> refPlane1;
        vec3<double> refPlane2;
        int boundStateIndex;
        std::tie(plane1, plane2, boundStateIndex) = calculatePlanes(part1, part2,
                                                                    const_cast<vec3<double> &>(pos1virtual));
        if (boundStateIndex > -1 ) {
            std::tie(refPlane1, refPlane2) = refPlanes[boundStateIndex];
            // Torque to align plane1:
            double cosSquared = (plane1 * refPlane1 + 1) * (plane1 * refPlane1 + 1);
            double cosThird = cosSquared * (plane1 * refPlane1 + 1);
            auto derivativeAngularPotential = (angularStrength / 4.0) * cosThird * plane1.cross(refPlane1);
            torque1 += derivativeAngularPotential; // Plus sign since plane1 x plane2 defined torque in particle 1
            torque2 -= derivativeAngularPotential;
            // Torque to align plane2:
            cosSquared = (plane2 * refPlane2 + 1) * (plane2 * refPlane2 + 1);
            cosThird = cosSquared * (plane2 * refPlane2 + 1);
            derivativeAngularPotential = (angularStrength / 4.0) * cosThird * plane2.cross(refPlane2);
            torque1 += derivativeAngularPotential; // Plus sign since plane1 x plane2 defined torque in particle 1
            torque2 -= derivativeAngularPotential;
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

    /* Given two particles, returns planes (unit vectors) to be aligned against reference planes
     * calculates in setBoundStates function. These are very unique to the specific application. It also
     * returns and integer corresponding to the most likely bound state (from which the planes
     * are chosen). Note returns zero if patches are not close enough. */
    std::tuple<vec3<double>, vec3<double>, int> patchyParticleSTMV::calculatePlanes(particle &part1,particle &part2,
            vec3<double> &pos1virtual) {
        vec3<double> plane1;
        vec3<double> plane2;
        std::vector<vec3<double>> part1Patches;
        std::vector<vec3<double>> part2Patches;

        for (auto &patch : patchesCoordinates) {
            part2Patches.push_back(msmrdtools::rotateVec(patch, part2.orientation));
        }

        // Calculate positions of patches (we can assume the frame of ref fixed at origin in base orientation)
        std::vector<vec3<double>> part1PatchPositions = patchesCoordinates;
        std::vector<vec3<double>> part2PatchPositions;
        for(auto &patch : part2Patches) {
            auto relpos = part2.position - pos1virtual;
            auto originPatch = relpos + msmrdtools::rotateVec(patch, part1.orientation.conj());
            part2PatchPositions.push_back(originPatch);
        }

        /* Calculate all the relevant distances between patches */
        std::vector<double> patchesDistances;
        for (auto part1Patch : part1PatchPositions) {
            for (auto part2Patch : part2PatchPositions) {
                patchesDistances.push_back((part1Patch - part2Patch).norm());
            }
        }

        /* Find minimum distance of all possible patches pairs */
        auto minIndex = static_cast<int> (std::min_element(patchesDistances.begin(), patchesDistances.end()) -
                                          patchesDistances.begin());

        /* Obtain indexes of closest interacting patches*/
        auto part1PatchIndex = std::floor(minIndex/5);
        int part2PatchIndex = minIndex%5;
        int boundStateIndex = -1;

        // Only allow for orientational potential if the molecules are bound or almost bound
        if (patchesDistances[minIndex] >= angularPotentialTolerance) {
            plane1 = vec3<double>(0,0,0);
            return std::make_tuple(plane1, plane1, -1);
        }

        /* Check the closest patches match one of the interactions (i1:i2, its inverse, i3:i4,
         * its inverse and i5:i5) */
        if (part1PatchIndex == 0 and part2PatchIndex == 1) {
            boundStateIndex = 0;
        } else if (part1PatchIndex == 1 and part2PatchIndex == 0) {
            boundStateIndex = 1;
        } else if (part1PatchIndex == 2 and part2PatchIndex == 3) {
            boundStateIndex = 2;
        } else if (part1PatchIndex == 3 and part2PatchIndex == 2) {
            boundStateIndex = 3;
        } else if (part1PatchIndex == 4 and part2PatchIndex == 4) {
            boundStateIndex = 4;
        }

        // If interaction is not allowed return (0 to produce zero torque)
        if (boundStateIndex == -1 ) {
            plane1 = vec3<double>(0,0,0);
            return std::make_tuple(plane1, plane1, 0);
        }

        // Calculate planes to be aligned against reference planes
        plane1 = part2PatchPositions[part2PatchIndex] - part1PatchPositions[part1PatchIndex];
        plane2 = plane1.cross(part2PatchPositions[part2PatchIndex]);
        plane1 = plane1/plane1.norm();
        plane2 = plane2/plane2.norm();

        return std::make_tuple(plane1, plane2, boundStateIndex);
    }

    // Returns actual position of different parts of particle (m1 to m3 or i1 to i5). Used for pybind.
    std::array<double, 3> patchyParticleSTMV::getPartPosition(std::string particlePart, particle &part) {
        vec3<double> position = vec3<double>();
        std::array<double, 3> outputPosition;
        if (particlePart == "m1"){
            position = part.position + msmrdtools::rotateVec(structureCoordinates[0], part.orientation);
        }
        else if (particlePart == "m2"){
            position = part.position + msmrdtools::rotateVec(structureCoordinates[1], part.orientation);
        }
        else if (particlePart == "m3"){
            position = part.position + msmrdtools::rotateVec(structureCoordinates[2], part.orientation);
        }
        else if (particlePart == "i1"){
            position = part.position + msmrdtools::rotateVec(patchesCoordinates[0], part.orientation);
        }
        else if (particlePart == "i2"){
            position = part.position + msmrdtools::rotateVec(patchesCoordinates[1], part.orientation);
        }
        else if (particlePart == "i3"){
            position = part.position + msmrdtools::rotateVec(patchesCoordinates[2], part.orientation);
        }
        else if (particlePart == "i4"){
            position = part.position + msmrdtools::rotateVec(patchesCoordinates[3], part.orientation);
        }
        else if (particlePart == "i5"){
            position = part.position + msmrdtools::rotateVec(patchesCoordinates[4], part.orientation);
        }
        return {position[0],position[1], position[2]};
    };

    // Sets particle diameters (sigma). In case defaults values need to be overriden.
    void patchyParticleSTMV::setParticlesDiameters(double sigmaM1, double sigmaM2, double sigmaM3, double sigmaOfPatches) {
        sigma[0] = sigmaM1;
        sigma[1] = sigmaM2;
        sigma[2] = sigmaM3;
        sigmaPatches = sigmaOfPatches;
    };

    // Sets parameters for repulsive part of potential. In case defaults values need to be overriden.
    void patchyParticleSTMV::setRepulsivePotentialParameters(double epsilon, double a, double rstar) {
        epsRepulsive = epsilon;
        aRepulsive = a;
        rstarRepulsive = rstar;
    }

    // Sets parameters for interaction between patches part of potential. In case defaults values need to be overriden.
    void patchyParticleSTMV::setInteractingPatchesPotentialParameters(double epsilon, double a, double rstar,
                                                                 int interactionIndex) {
        epsPatchesList[interactionIndex] = epsilon;
        aPatchesList[interactionIndex] = a;
        rstarPatchesList[interactionIndex] = rstar;
    }




}
