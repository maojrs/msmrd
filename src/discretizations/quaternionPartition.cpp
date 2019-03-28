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
        halfSphericalPartition = std::make_unique<halfSpherePartition>(numSphericalSections);
    };

    std::tuple<std::vector<int>, std::vector<double> > quaternionPartition::makePartition() {
        double dtheta = 2*M_PI/numRadialSections; // from the vector part of the quaternion u*sin(theta/2)
        std::vector<double> radialSections;
        double drQuaternion;
        for (int i=0; i<numRadialSections; i++){
            // Transform to quaternion scaling
            drQuaternion = 2*std::asin(i*dtheta); // rescaling with vector part of quaternion u*sin(theta/2)
            radialSections.push_back(drQuaternion);
        }




    };


}
