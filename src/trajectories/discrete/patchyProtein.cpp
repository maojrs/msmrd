//
// Created by maojrs on 4/8/19.
//

#include "trajectories/discrete/patchyProtein.hpp"

namespace msmrd {
    /*
     * Chooses how to discretize the full trajectory (trajectoryData) into a discretized trajectory
     * to be analyzed and extracted into a Markov state model. This specific discretization will follow
     * the core MSM approach
     *
     * @param reltiveDistanceCutOff determines the cutoff for the relative position, after this value is exceeded,
     * the partition is no longer effective and the getSectionNumber function returns 0 (unbound state)
     * @param numSphericalSectionsPos number of sections in the surface partition of the unit sphere to discretize
     * the direction of the relative position.
     * @param numRadialSectionsQuat number of radial sections in volumetric 3D sphere partition used to
     * discretize the relative orientation quaternion
     * @param numSphericalSectionsQuat number of spherical sections in each shell of the quaternion discretization.
     * Each radial shell in the quaternion partition has numSphericalShells section in it.
     */
    patchyProtein::patchyProtein(unsigned long Nparticles, int bufferSize, double relativeDistanceCutOff,
                                 int numSphericalSectionsPos, int numRadialSectionsQuat,
                                 int numSphericalSectionsQuat) : trajectoryPositionOrientation(Nparticles, bufferSize) {
        discreteTrajectoryData.reserve(bufferSize);
        positionOrientationPart = std::make_unique<positionOrientationPartition>(relativeDistanceCutOff,
                numSphericalSectionsPos, numRadialSectionsQuat, numSphericalSectionsQuat);
        setMetastableRegions();
    }


    void patchyProtein::sampleDiscreteTrajectory(double time, std::vector<particle> &particleList) {

    }

    /*
     * Define metastable relative orientations (including symmetric quaternion)
     */
    void patchyProtein::setMetastableRegions() {
        // Define 2 main relative rotations (quaternions) that define bounds states A and B.
        const int numStates = 6;
//        std::array<vec3<double>, numStates> axisAngleVecs;
//        rotMetastableStates.resize(numStates);
//        axisAngleVecs[0] = vec3<double>(0, 0, M_PI - 3.0*M_PI/5);
//        axisAngleVecs[1] = vec3<double>(0, 0, M_PI);
//        for (int i = 0; i < numStates; i++){
//            rotMetastableStates[i] = msmrdtools::axisangle2quaternion(axisAngleVecs[i]);
//        }
    }


}
