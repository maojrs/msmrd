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
               angularStrength(angularStrength) {

        setPotentialParameters();

        //setBoundStates();
    };

    patchyProteinMarkovSwitch::patchyProteinMarkovSwitch(double sigma, double strength, double angularStrength,
                                 std::vector<std::vector<double>> patchesCoordinatesA,
                                 std::vector<std::vector<double>> patchesCoordinatesB)
            :  patchyProtein(sigma, strength, patchesCoordinatesA, patchesCoordinatesB),
               angularStrength(angularStrength) {

        setPotentialParameters();

        //setBoundStates();
    };




    /* Set potentials parameters for overall potential. Consists of isotropic attractive and
     * repulsive parts plus two types of patches interactions*/
    void patchyProteinMarkovSwitch::setPotentialParameters() {
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
        epsAttractive = -0.15*strength; //-0.05*strength;
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

        // Sets this specific parameter unique for this potential, see enableDisableMSM function.
        minimumR = 1.25;
    }

    /* Checks is MSM must be deactivated in certain particles. In this case, if bounded or close to bounded
     * disable MSM in particle withs type 1 and state 0. This needs to be hardcoded here for each example.
     * Not used at the moment. */
    void patchyProteinMarkovSwitch::enableDisableMSM(vec3<double>relPosition, particle &part1, particle &part2) {
        if (relPosition.norm() <= minimumR && part2.state == 0) {
            part2.deactivateResetMSM();
            part2.activeMSM = false;
        } else {
            part2.activateMSM();
        }
    }


    // Evaluates potential at given positions and orientations of two particles
    double patchyProteinMarkovSwitch::evaluate(particle &part1, particle &part2) {
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

        // Evaluate patches potential if close enough and if particle 2 is in state 0
        if (rvec.norm() <= 2*sigma and part2.state == 0 and patchesActive) {
            // Use default patches auxiliary parent function without angular dependence
            patchesPotential = evaluatePatchesPotential(part1, part2, pos1virtual, patchesCoords1, patchesCoords2);
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
    std::array<vec3<double>, 4> patchyProteinMarkovSwitch::forceTorque(particle &part1, particle &part2)  {

        // Calculate relative position
        std::array<vec3<double>, 2> relPos = relativePositionComplete(part1.position, part2.position);
        vec3<double> pos1virtual = relPos[0]; // virtual part1.position if periodic boundary; otherwise part1.position.
        vec3<double> rvec = relPos[1]; //part2.position - part1.position;

        // Calculate relative orientation
        auto relOrientation = part2.orientation * part1.orientation.conj();

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
         * disable and reset the MSM on particle2 when it is on state 0. This is not used for the MSM/RD example, but
         * it is left for possible future implementations. */
        //enableDisableMSM(rvec, part1, part2);

        // Check if particle is in any of the bound states. Not in use, left here for possible future implementations
        //auto particlesBound = isBound(rvec, relOrientation);
        // If bound, deactivate the unboundMSM, otherwise, keep it active (assumes unbound MSM only of particle 2.).
        //if (particlesBound) {
        //    part2.deactivateMSM();
        //} else {
        //    part2.activateMSM();
        //}

        // Calculate forces and torque due to patches interaction, if close enough and if particle 2 is in state 0
        if ( rvec.norm() <= 2*sigma and part2.state == 0 and patchesActive) {
            // Calculate forces and torque due to patches interaction using auxiliary function
            auto forcTorqPatches = forceTorquePatches(part1, part2, pos1virtual, patchesCoords1, patchesCoords2);
            auto force1 = forcTorqPatches[0];
            auto torque1 = forcTorqPatches[1];
            auto force2 = forcTorqPatches[2];
            auto torque2 = forcTorqPatches[3];

            /* Explicit angular dependence. */
            // Get planes needed to be aligned by torque
            vec3<double> plane1;
            vec3<double> plane2;
            std::tie(plane1, plane2) = calculatePlanes(part1, part2, patchesCoords1, patchesCoords2);
            /* If one uses only one minima, use potential of -(1/256)*(cos(theta) + 1)^8, with theta the angle
             * between the unitary vector of each plane */
            double cosSquared = (plane1 * plane2 + 1) * (plane1 * plane2 + 1);
            double cosSeventh = cosSquared * cosSquared * cosSquared * (plane1 * plane2 + 1);
            auto derivativeAngluarPotential = 0.5*angularStrength*sigma*(8.0/256.0) * cosSeventh * plane1.cross(plane2);
            torque1 += derivativeAngluarPotential; // Plus sign since plane1 x plane2 defined torque in particle 1
            torque2 -= derivativeAngluarPotential;

            return {force + force1, torque1, -1.0*force + force2, torque2};
        } else {
            auto zeroTorque = vec3<double> (0.0, 0.0, 0.0);
            return {force, zeroTorque, -1.0 * force, zeroTorque};
        }
    }


    /* Given two quaternions/orientations, returns planes(unit vectors) to be aligned by torque. These
     * may vary depending on the physical arrangement of your molecules and how the patches are ordered and
     * selected. Currently set up for particle 1 with 6 binding patches and particle 2 only with one.*/
    std::tuple<vec3<double>, vec3<double>> patchyProteinMarkovSwitch::calculatePlanes(particle &part1,
                                                                 particle &part2,
                                                                 const std::vector<vec3<double>> patches1,
                                                                 const std::vector<vec3<double>> patches2) {

        std::vector<vec3<double>> part1PatchNormals;
        std::vector<vec3<double>> part2PatchNormals;
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

        /* Find minimum distance and match corresponding expected orientation (0, 1, 2,...5 in original implementation)
         * depending on the patches that are closer to each other) */
        auto minIndex = static_cast<int> (std::min_element(patchesDistances.begin(), patchesDistances.end()) -
                                          patchesDistances.begin());

        /* Assign using indexDict second patch in patchesCoordinatesA to do cross product with
         * patchesCoordinatesA[minindex]. Dependent on implemenatiation. In this case, assume the patches
         * follow this order: patchesCoordinatesA = [1.,0.,0.], [0.,1.,0.], [0.,0.,1.],
         * [-1.,0.,0.], [0.,-1.,0.], [0.,0.,-1.]. Note the orientations should match the discrete trajectory
         * definition; in this case, the one given by setBoundStates() in trajectories/discrete/patchyProtein.*/
        std::array<int, 6> indexDict = {1, 3, 4, 4, 0, 4};
        int secondIndex = indexDict[minIndex];

        /* It further assumes, patchesCoordinatesB = [1.,0.,0.]; defining patchRef (orthogonal to previous patch), they
         * define the orientation of particle 2.  */
        vec3<double> patchRef = vec3<double> (0.0, 1.0, 0.0);

        // Calculate unitary vectors describing planes where particle center and first two patches are
        vec3<double> plane1;
        vec3<double> plane2;
        auto patchRefRotated = msmrdtools::rotateVec(patchRef, part2.orientation);
        plane1 = part1PatchNormals[minIndex].cross(part1PatchNormals[secondIndex]); // should always be perpendicular.
        plane2 = part2PatchNormals[0].cross(patchRefRotated); // should always be perpendicular.

        // Not neccesary to renormalize; they must always be normal (or close enough up to machine-eps)
        //plane1 = plane1/plane1.norm();
        //plane2 = plane2/plane2.norm();

        return std::make_tuple(plane1, plane2);
    }


    /* Copied from discreteTrajectory to identify bound states. Checks if particle is in any of the patchyProtein
     * trajectory bound states. Not in use at the moment */
    bool patchyProteinMarkovSwitch::isBound(vec3<double> relativePosition, quaternion<double> relativeOrientation) {
        // Check if it matches a bound states, if so return true otherwise return false.
        vec3<double> relPosCenter;
        quaternion<double> relQuatCenter;
        for (int i = 0; i < 6; i++) {
            relPosCenter = std::get<0>(boundStates[i]);
            relQuatCenter = std::get<1>(boundStates[i]);

            if ( (relPosCenter - relativePosition).norm() <= 0.12) {
                auto angleDistance = msmrdtools::quaternionAngleDistance(relQuatCenter, relativeOrientation);
                if  (angleDistance < 0.12*2*M_PI) {
                    return true;
                }
            }
        }
        return false;
    }


    /* Copied from patchyProteinTrajectory. Sets bound states (metastable regions) of the patchy protein
     * implementation. Should match those of the patchyProtein trajectory. The centers of the
     * metastable regions are given by a tuple of relative position and relative orientation. The size of the
     * regions are determined by tolerancePosition and toleranceOrientation. Not required at the moment. */
    void patchyProteinMarkovSwitch::setBoundStates() {
        /* Define relative position vectors from particle 1 at the origin. These two patches
         * point in the same direction as the two patches in the dimer. */
        std::array<vec3<double>, 6> relPos;
        relPos[0] = {1., 0., 0.};
        relPos[1] = {0., 1., 0.};
        relPos[2] = {0., 0., 1.};
        relPos[3] = {-1., 0., 0.};
        relPos[4] = {0., -1., 0.};
        relPos[5] = {0., 0., -1.};
        /* Relative rotations (assuming particle 1 fixed) of particle 2 that yield the 6 bound states
         * in the axis-angle representation. (One needs to make drawing to understand)*/
        std::array<vec3<double>, 6> rotations;
        rotations[0] = {0.0, 0.0, M_PI}; //ok
        rotations[1] = {0.0, 0.0, -M_PI / 2.0}; //ok
        rotations[2] = {0.0, M_PI / 2.0, 0.0}; //ok
        rotations[3] = {0.0, 0.0, 0.0}; //ok
        rotations[4] = {0.0, 0.0, M_PI / 2.0}; //ok
        rotations[5] = {0.0, -M_PI / 2.0, 0.0}; //ok
        /*Convert rotations in the axis angle representation to quaternions */
        std::array<quaternion<double>, 6> quatRotations;
        for (int i = 0; i < 6; i++) {
            quatRotations[i] = msmrdtools::axisangle2quaternion(rotations[i]);
        }
        /* Fill bound states with corresponding combinations of relative position vectors and quaternion orientations.
         * Note we need to take the conj, so it matches the relative orientation from particle 1, as used to define
         * the states in trajectories/discrete/discreteTrajectory.hpp */
        boundStates[0] = std::make_tuple(relPos[0], quatRotations[0].conj());
        boundStates[1] = std::make_tuple(relPos[1], quatRotations[1].conj());
        boundStates[2] = std::make_tuple(relPos[2], quatRotations[2].conj());
        boundStates[3] = std::make_tuple(relPos[3], quatRotations[3].conj());
        boundStates[4] = std::make_tuple(relPos[4], quatRotations[4].conj());
        boundStates[5] = std::make_tuple(relPos[5], quatRotations[5].conj());
    }


}

