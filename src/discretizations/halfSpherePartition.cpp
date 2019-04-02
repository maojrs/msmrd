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
        scaling = 2;
        partitionSphere();
    }

}