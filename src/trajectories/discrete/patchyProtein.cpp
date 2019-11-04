//
// Created by maojrs on 4/8/19.
//

#include "trajectories/discrete/patchyProtein.hpp"

namespace msmrd {

    /*
     * Constructors and definition of bound states of patchy protein classes.
     * IMPORTANT: Need to define pointer to positionOrientationPartition and
     * call setBoundStates in constructors of childs of discreteTrajectory.
     */

    patchyProtein::patchyProtein(unsigned long Nparticles, int bufferSize) :
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

    patchyProtein::patchyProtein(unsigned long Nparticles, int bufferSize, double rLowerBound, double rUpperBound) :
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
    void patchyProtein::setBoundStates() {
        double angleDiff = 3 * M_PI / 5; // angle difference to form a pentamer
        /* Define relative position vectors from particle 1 at the origin. These two patches
         * point in the same direction as the two patches in the dimer. */
        vec3<double> relPos1 = {1., 0., 0.};
        vec3<double> relPos2 = {0., 1., 0.};
        vec3<double> relPos3 = {0., 0., 1.};
        vec3<double> relPos4 = {-1., 0., 0.};
        vec3<double> relPos5 = {0., -1., 0.};
        vec3<double> relPos6 = {0., 0., -1.};
        /* Relative rotations (from particle 1) of particle 2 that yield the 6 bound states
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
        // Fill bound states with corresponding combinations of relative position vectors and quaternion orientations.
        boundStates[0] = std::make_tuple(relPos1, quatRotations[0]);
        boundStates[1] = std::make_tuple(relPos2, quatRotations[1]);
        boundStates[2] = std::make_tuple(relPos3, quatRotations[2]);
        boundStates[3] = std::make_tuple(relPos4, quatRotations[3]);
        boundStates[4] = std::make_tuple(relPos5, quatRotations[4]);
        boundStates[5] = std::make_tuple(relPos6, quatRotations[5]);
    }


}
