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
        radialSections.resize(numRadialSections+1);
        numTotalSections = numSphericalSections*(numRadialSections -1) + 1;
        makeRadialPartition();
    };

    // Defines the location of the radial cuts between the origin and r=1.
    void quaternionPartition::makeRadialPartition() {
        double dr = 1.0/numRadialSections;
        radialSections[0] = 0.0;
        for (int i=1; i<=numRadialSections; i++){
            radialSections[i] = radialSections[i-1] + dr;
        }
    };

    // Obtains the section number in volumetric half sphere partition, given a coordinate within the half unit sphere.
    int quaternionPartition::getSectionNumber(quaternion<double> quatCoordinate) {
        int sectionNumber;
        double r = quatCoordinate.norm();
        if (r > 1.0005) {
            throw std::range_error("Unit quaternion cannot be larger than one");
        }
        // Reduce quaternion coordinate from (s,x,y,z) to (s,x,z)
        vec3<double> reducedCoordinate{quatCoordinate[0], quatCoordinate[1], quatCoordinate[3]};
        double rReduced = reducedCoordinate.norm();

        for (int i = 0; i < numRadialSections; i++){
            if (rReduced <= radialSections[i+1]){
                if (i == 0) {
                    sectionNumber = 1;
                    break;
                } else {
                    sectionNumber = numSphericalSections*(i-1) + 1;
                    sectionNumber += sphericalPartition->getSectionNumber(reducedCoordinate);
                    break;
                }

            }
        }
        return sectionNumber;
    };

    /* Gets volumetric interval delimiter of the section corresponding to secNumber. The function return three
     * intervals, one in the r direction, one in the phi angle (polar) and one in the theta angle (azimuthal). */
    std::tuple<std::vector<double>, std::vector<double>,
            std::vector<double>> quaternionPartition::getSectionIntervals(int secNumber) {
        std::vector<double> rInterval(2, 0);
        std::vector<double> phiInterval(2,0);
        std::vector<double> thetaInterval(2,0);
        // If on first section, return the full inner sphere interval limits.
        if (secNumber == 1) {
            rInterval = {0.0, radialSections[1]};
            phiInterval = {0.0, M_PI};
            thetaInterval = {0, 2*M_PI};

        } else {
            /* Otherwise, find on which shell is the particle and use the halSphericalPartition to obtain
             * the correct angular intervals */
            std::tuple<std::vector<double>, std::vector<double>> angles;
            double rindex = -1;
            int sectionNumModule = (secNumber - 1) % numSphericalSections;
            if (sectionNumModule == 0) {
                rindex = static_cast<int> (std::floor((secNumber - 1) / numSphericalSections));
                angles = sphericalPartition->getAngles(numSphericalSections);
            } else {
                rindex = static_cast<int> (std::floor((secNumber - 1) / numSphericalSections) + 1);
                angles = sphericalPartition->getAngles(sectionNumModule);
            }
            rInterval = {radialSections[rindex], radialSections[rindex + 1]};
            phiInterval = std::get<0>(angles);
            thetaInterval = std::get<1>(angles);
        }
        return std::make_tuple(rInterval, phiInterval, thetaInterval);
    };

}
