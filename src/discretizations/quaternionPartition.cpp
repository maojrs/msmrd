//
// Created by maojrs on 3/27/19.
//

#include "discretizations/quaternionPartition.hpp"

namespace msmrd {

    /*
     * Constructor creates a spherical equal area partition on the surface of a sphere.
     * @param numSections the number of sections that the sphere should be partioned in.
     */
    quaternionPartition::quaternionPartition(int numRadialSections, int numSphericalSections):
            numRadialSections(numRadialSections), numSphericalSections(numSphericalSections){
        sphericalPartition = std::make_unique<spherePartition>(numSphericalSections);
    };

    std::tuple<std::vector<int>, std::vector<double> > quaternionPartition::makePartition() {
        double dtheta = 2*M_PI/numRadialSections;
        std::vector<double> radialSections;

        for (int i=0; i<numRadialSections; i++){
            // Transform to quaternion scaling
            double drQuaternion = 2*std::asin(i*dtheta);
            radialSections.push_back(drQuaternion);
        }




    };


}
