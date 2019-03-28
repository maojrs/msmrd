//
// Created by maojrs on 3/28/19.
//
#include "discretizations/halfSpherePartition.hpp"

namespace msmrd {

    /*
     * Constructor creates a spherical equal area partition on the surface of a sphere.
     * @param numSections the number of sections that the sphere should be partioned in.
     */
    halfSpherePartition::halfSpherePartition(int numSecs) {
        numSections = numSecs;
        double scaling = 2; // Scales to half a sphere, other scaling must be carefully handled
        auto partition = partitionSphere(numSections, scaling);
        regionsPerCollar = std::get<0>(partition);
        phis = std::get<1>(partition);
        thetas = std::get<2>(partition);
    }

}