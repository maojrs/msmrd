//
// Created by maojrs on 1/22/19.
//
#include "trajectories/discrete/patchyDimer.hpp"

namespace msmrd {

    /*
     * Constructors and definition of bound states of patchy dimer classes.
     * IMPORTANT: Need to define pointer to positionOrientationPartition and
     * call setBoundStates in constructors of childs of discreteTrajectory.
     */

    /*
     * Patchy dimer (the original classs). The first toy model to test MSM/RD. Consists of 8 bound states,
     * corresponding to all the possible ways to bound through the two patches of each particle combined
     * with two stable angular configurations. Use with patchyParticleAngular potential with numAngularMinima = 2
     */

    patchyDimer::patchyDimer(unsigned long Nparticles, int bufferSize) :
    discreteTrajectory(Nparticles, bufferSize) {

        setRadialBounds(1.25, 2.25);
        setTolerances(0.12, 0.12*2*M_PI);

        int numSphericalSectionsPos = 7; //7; //7;
        int numRadialSectionsQuat = 5; // 3; //5;
        int numSphericalSectionsQuat =7; // 6; //7;
        positionOrientationPart = std::make_unique<positionOrientationPartition>(rUpperBound, numSphericalSectionsPos,
                                                                                 numRadialSectionsQuat,
                                                                                 numSphericalSectionsQuat);
        setBoundStates();
    };

    patchyDimer::patchyDimer(unsigned long Nparticles, int bufferSize, double rLowerBound, double rUpperBound) :
        discreteTrajectory(Nparticles, bufferSize, rLowerBound, rUpperBound) {

        setRadialBounds(rLowerBound, rUpperBound);
        setTolerances(0.12, 0.12*2*M_PI);

        int numSphericalSectionsPos = 7; //7; //7;
        int numRadialSectionsQuat = 5; // 3; //5;
        int numSphericalSectionsQuat =7; // 6; //7;
        positionOrientationPart = std::make_unique<positionOrientationPartition>(rUpperBound, numSphericalSectionsPos,
                                                                                 numRadialSectionsQuat,
                                                                                 numSphericalSectionsQuat);
        setBoundStates();
    };



    /* Sets bound states (metastable regions) of this patchy dimer implementation (two equal patcy particles,
     * each with two patches an angle angleDiff away with two stable relative orientations). The centers of the
     * metastable regions are given by a tuple of relative position and relative orientation. The size of the
     * regions are determined by tolerancePosition and toleranceOrientation*/
    void patchyDimer::setBoundStates() {
        double angleDiff = 3 * M_PI / 5; // angle difference to form a pentamer
        /* Define relative position vectors from particle 1 at the origin. These two patches
         * point in the same direction as the two patches in the dimer. */
        vec3<double> relPos1 = {std::cos(angleDiff / 2.0), std::sin(angleDiff / 2.0), 0};
        vec3<double> relPos2 = {std::cos(angleDiff / 2.0), std::sin(-angleDiff / 2.0), 0};
        vec3<double> relPos1orthogonal = {-1.0 * std::sin(angleDiff / 2.0), std::cos(angleDiff / 2.0), 0.0};
        vec3<double> relPos2orthogonal = {std::sin(angleDiff / 2.0), std::cos(angleDiff / 2.0), 0.0};
        /* Relative rotations (from particle 1) of particle 2 that yield the 8 bound states
         * in the axis-angle representation. (One needs to make drawing to understand)*/
        std::array<vec3<double>, 8> rotations;
        rotations[0] = M_PI * relPos1orthogonal; //ok
        rotations[1] = {0.0, 0.0, -2 * M_PI / 5.0}; //ok
        rotations[2] = {0.0, 0.0, M_PI}; //ok
        rotations[3] = {0.0, M_PI, 0.0}; //ok
        // --first 4 rotations correspond to binding on top patch of particle 1, next 4 rotations to bottom patch
        rotations[4] = M_PI * relPos2orthogonal; //ok
        rotations[5] = {0.0, 0.0, 2 * M_PI / 5.0}; //ok
        rotations[6] = {0.0, 0.0, M_PI}; //ok
        rotations[7] = {0.0, M_PI, 0.0}; //ok
        /*Convert rotations in the axis angle representation to quaternions */
        std::array<quaternion<double>, 8> quatRotations;
        for (int i = 0; i < 8; i++) {
            quatRotations[i] = msmrdtools::axisangle2quaternion(rotations[i]);
        }
        // Fill bound states with corresponding combinations of relative position vectors and quaternion orientations.
        boundStates[0] = std::make_tuple(relPos1, quatRotations[0]);
        boundStates[1] = std::make_tuple(relPos1, quatRotations[1]);
        boundStates[2] = std::make_tuple(relPos1, quatRotations[2]);
        boundStates[3] = std::make_tuple(relPos1, quatRotations[3]);
        boundStates[4] = std::make_tuple(relPos2, quatRotations[4]);
        boundStates[5] = std::make_tuple(relPos2, quatRotations[5]);
        boundStates[6] = std::make_tuple(relPos2, quatRotations[6]);
        boundStates[7] = std::make_tuple(relPos2, quatRotations[7]);
    }



    /*
     * Patchy dimer 2. Base model to obtain MSM for MSM/RD of pentamer formation. Unlike patchy dimer class, this
     * one only has one stable angular configuration, yielding only 4 bound states. Use with patchyParticleAngular
     * potential with numAngularMinima = 2
     */

    patchyDimer2::patchyDimer2(unsigned long Nparticles, int bufferSize) :
            discreteTrajectory(Nparticles, bufferSize) {

        setRadialBounds(1.25, 2.25);
        setTolerances(0.12, 0.12*2*M_PI);

        int numSphericalSectionsPos = 7; //7; //7;
        int numRadialSectionsQuat = 5; // 3; //5;
        int numSphericalSectionsQuat =7; // 6; //7;
        positionOrientationPart = std::make_unique<positionOrientationPartition>(rUpperBound, numSphericalSectionsPos,
                                                                                 numRadialSectionsQuat,
                                                                                 numSphericalSectionsQuat);
        setBoundStates();

    };

    patchyDimer2::patchyDimer2(unsigned long Nparticles, int bufferSize, double rLowerBound, double rUpperBound) :
            discreteTrajectory(Nparticles, bufferSize, rLowerBound, rUpperBound) {

        setRadialBounds(rLowerBound, rUpperBound);
        setTolerances(0.12, 0.12*2*M_PI);

        int numSphericalSectionsPos = 7; //7; //7;
        int numRadialSectionsQuat = 5; // 3; //5;
        int numSphericalSectionsQuat =7; // 6; //7;
        positionOrientationPart = std::make_unique<positionOrientationPartition>(rUpperBound, numSphericalSectionsPos,
                                                                                 numRadialSectionsQuat,
                                                                                 numSphericalSectionsQuat);
        setBoundStates();
    };


    /* Sets bound states (metastable regions) of this patchy dimer implementation (two equal patcy particles,
     * each with two patches an angle angleDiff away and with one stable relative orientation). The centers of
     * the metastable regions are given by a tuple of relative position and relative orientation. The size of
     * the regions are determined by tolerancePosition and toleranceOrientation*/
    void patchyDimer2::setBoundStates() {
        double angleDiff = 3 * M_PI / 5; // angle difference to form a pentamer
        /* Define relative position vectors from particle 1 at the origin. These two patches
         * point in the same direction as the two patches in the dimer. */
        vec3<double> relPos1 = {std::cos(angleDiff / 2.0), std::sin(angleDiff / 2.0), 0};
        vec3<double> relPos2 = {std::cos(angleDiff / 2.0), std::sin(-angleDiff / 2.0), 0};
        vec3<double> relPos1orthogonal = {-1.0 * std::sin(angleDiff / 2.0), std::cos(angleDiff / 2.0), 0.0};
        vec3<double> relPos2orthogonal = {std::sin(angleDiff / 2.0), std::cos(angleDiff / 2.0), 0.0};
        /* Relative rotations (from particle 1) of particle 2 that yield the 4 bound states
         * in the axis-angle representation. (One needs to make drawing to understand)*/
        std::array<vec3<double>, 4> rotations;
        rotations[0] = M_PI * relPos1orthogonal; // part1Patch1 with part2patch1
        rotations[1] = {0.0, 0.0, -2 * M_PI / 5.0}; // part1Patch1 with part2patch2
        // --first 2 rotations correspond to binding on top patch (1) of part1, next 2 rotations to bottom patch (2).
        rotations[2] = {0.0, 0.0, 2 * M_PI / 5.0}; // part1Patch2 with part2patch1
        rotations[3] = M_PI * relPos2orthogonal; // part1Patch2 with part2patch2
        /*Convert rotations in the axis angle representation to quaternions */
        std::array<quaternion<double>, 4> quatRotations;
        for (int i = 0; i < 4; i++) {
            quatRotations[i] = msmrdtools::axisangle2quaternion(rotations[i]);
        }
        // Fill bound states with corresponding combinations of relative position vectors and quaternion orientations.
        boundStates[0] = std::make_tuple(relPos1, quatRotations[0]);
        boundStates[1] = std::make_tuple(relPos1, quatRotations[1]);
        boundStates[2] = std::make_tuple(relPos2, quatRotations[2]);
        boundStates[3] = std::make_tuple(relPos2, quatRotations[3]);
    }

}