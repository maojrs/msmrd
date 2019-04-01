#pragma once
#include <vector>
#include <tuple>
#include <memory>
#include "quaternion.hpp"
#include "spherePartition.hpp"
#include "vec3.hpp"


namespace msmrd {
    /*
     * This class creates a partition on the volume of a unit sphere. Its goal is to partition quaternion space
     * (rotations), which corresponds to a unit four-dimensional half sphere. This space can be projected into a 3D
     * volume inside the unit sphere for easier discretization (similar to a stereographic projection on the top half
     * 3D sphere). It uses spherePartition.hpp to divide the 3d sphere. In order to obtain the half unit quaternion
     * sphere projection of (s,ix,jy,kz), one of the imaginary variables must be flattened (x, y or z); e.g. the new
     * coordinates could simply be (s,x,z) and recover y with y=+sqrt(1-ss-xx-zz). If the s variable is flatten,
     * important information about the rotation would be lost. On the other hand any quaternion with negative y
     * (x or z), can be rewritten as a quaternion with positive y (x or z), so no information of the rotation is
     * lost when cutting negative values of y (x or z). (Quaternion symmetries only along imaginary variables)
     *
     * The recommended quaternion coordinates to do the 3D discretization are (s,x,z), while the y coordinate is
     * flattened. Half sphere lives in y>=0, while negative y values should be rewritten as rotations with y>0 (always
     * possible). The y coordinate is recovered as y=+sqrt(1-s*s-x*x-z*z). Note the choice of y was arbitrary. it could
     * have been x or z as well.
     */
    class quaternionPartition {
    public:
        int numRadialSections;
        int numSphericalSections;
        int numTotalSections = 0;
        std::unique_ptr<spherePartition> sphericalPartition;
        std::vector<double> radialSections;
        /**
         * @param numRadialSections number of radial sections in volumetric 3D sphere partition
         * @param numSphericalSections number of spherical sections in surface of sphere. Each radial shell in
         * the partition has numSphericalShells section in it.
         * @param numTotalSections total number of sections in partition
         * @param sphericalPartition pointer to pherical partition on the surface. Along with the
         * radial sections it is used to define the volumetric sections.
         * @param radialSections vector with radial values to define cuts of radial sections. Along with the
         * halfSphericalPartition it is used to define the volumetric sections.
         */

        quaternionPartition(int numRadialSections, int numSphericalSections);

        // Defines the location of the radial cuts in partition.
        void makeRadialPartition();

        /* Gets the section number in the volumetric partition of the half sphere given a coordinate
         * inside the half unit sphere. Uses (s,x,z) as reduce coordinate */
        int getSectionNumber(quaternion<double> quatCoordinate);

        /* Gets volumetric interval delimiter of the section corresponding to secNumber. The function return three
        * intervals, one in the r direction, one in the phi angle (polar) and one in the theta angle (azimuthal). */
        std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> getSectionIntervals(int secNumber);
    };

}