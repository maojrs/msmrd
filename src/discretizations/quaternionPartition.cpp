//
// Created by maojrs on 3/27/19.
//

#include "discretizations/quaternionPartition.hpp"

namespace msmrd {

    /*
     * Constructor creates a volumetric partition of the unit sphere in 3D that discretizes a unit quaternion.
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

    /* Gets the section number in the volumetric partition of the half sphere given a coordinate
     * inside the half unit sphere. Uses (s,x,z) as reduce coordinate */
    int quaternionPartition::getSectionNumber(quaternion<double> quatCoordinate) {
        int sectionNumber;
        double r = quatCoordinate.norm();
        if (r > 1.0005) {
            throw std::range_error("Unit quaternion cannot be larger than one");
        }
        /* Reduce quaternion coordinate from (s,x,y,z) to (x,y,z), uses the fact
         * that rotation of q = rotation of -q */
        vec3<double> reducedCoordinate;
        if (quatCoordinate[0] >= 0) {
            reducedCoordinate = {quatCoordinate[1], quatCoordinate[2], quatCoordinate[3]};
        } else {
            reducedCoordinate = {-1.0*quatCoordinate[1], -1.0*quatCoordinate[2], -1.0*quatCoordinate[3]};
        }
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
    std::tuple<std::array<double, 2>, std::array<double, 2>,
            std::array<double, 2>> quaternionPartition::getSectionIntervals(int secNumber) {
        std::array<double, 2> rInterval;
        std::array<double, 2> phiInterval;
        std::array<double, 2> thetaInterval;
        // If on first section, return the full inner sphere interval limits.
        if (secNumber == 1) {
            rInterval = {0.0, radialSections[1]};
            phiInterval = {0.0, M_PI};
            thetaInterval = {0, 2*M_PI};

        } else {
            /* Otherwise, find on which shell is the particle and use the sphericalPartition to obtain
             * the correct angular intervals */
            std::tuple<std::array<double, 2>, std::array<double, 2>> angles;
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

    // Gets full partition: radialSections, regionsPerCollar, phis and thetas
    std::tuple<std::vector<double>, std::vector<int>,
            std::vector<double>, std::vector<std::vector<double>>> quaternionPartition::getPartition() {
        auto regionsPerCollar = sphericalPartition->regionsPerCollar;
        auto phis = sphericalPartition->phis;
        auto thetas = sphericalPartition->thetas;
        return std::make_tuple(radialSections, regionsPerCollar , phis, thetas);
    };

    // Pybind version of getSectionNumber, uses arbitrary vectors instead of vec3.
    int quaternionPartition::getSectionNumberPyBind(std::vector<double> coord) {
        return getSectionNumber(quaternion<double>{coord[0], coord[1], coord[2], coord[3]});
    }

    // Returns total number of sections, num radial sections and number of spherical sections
    std::vector<int> quaternionPartition::getNumSections(){
        return std::vector<int>{numTotalSections, numRadialSections, numSphericalSections};
    }

}
