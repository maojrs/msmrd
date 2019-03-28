#pragma once
#include <vector>
#include <tuple>
#include <memory>
#include "vec3.hpp"
#include "quaternion.hpp"
#include "halfSpherePartition.hpp"


namespace msmrd {
    /*
     * This class creates a partition on the volume of a half unit sphere. Its goal is to partition quaternion space,
     * which corresponds to a unit four-dimensional half sphere. This space can be projected into a 3D volume inside
     * the unit sphere for easier discretization (similar to a stereographic projection on the top half sphere). It
     * uses spherePartition.hpp to divide the 3d sphere.
     */
    class quaternionPartition {
    public:
        int numRadialSections;
        int numSphericalSections;
        std::unique_ptr<halfSpherePartition> halfSphericalPartition;

        quaternionPartition(int numRadialSections, int numSphericalSections);

        std::tuple<std::vector<int>, std::vector<double> > makePartition();
    };

}