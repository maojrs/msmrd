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
     * sphere projection of (s,ix,jy,kz), one of the variables must be flattened (s, x, y or z); e.g. the new
     * coordinates could simply be (x,y,z) and recover s with s=+sqrt(1-xx-yy-zz). This uses the fact that the
     * quaternion q represents the same rotation as quaternion -q.
     *
     * The recommended quaternion coordinates to do the 3D discretization are (x,y,z), while the s coordinate is
     * flattened. Half sphere lives in s>=0, while negative s values should be rewritten as rotations with s>0 (always
     * possible). Note the choice of s was arbitrary. it could have been x,y or z as well.
     *
     * Note section numbering (secNumber) starts in 1 and not zero.
     */
    class quaternionPartition {
    protected:

        // Defines the location of the radial cuts in partition.
        void makeRadialPartition();

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

        int getSectionNumber(quaternion<double> quatCoordinate);

        std::tuple<std::array<double, 2>, std::array<double, 2>,
                std::array<double, 2>> getSectionIntervals(int secNumber);



        /* Other not so important functions (mostly for PyBindings)*/

        std::tuple<std::vector<double>, std::vector<int>,
                std::vector<double>, std::vector<std::vector<double>>> getPartition();

        int getSectionNumberPyBind(std::vector<double> coord);

        std::vector<int> getNumSections();




    };

}