//
// Created by maojrs on 3/28/19.
//
#pragma once
#include "spherePartition.hpp"

namespace msmrd {
    /*
     * This class creates an equal area partition on the surface of a sphere. It is a c++ copy and extension
     * of the python code in module msmrd2.tools.spherePartition.
     *
     * Note section numbering (secNumber) starts in 1 and not zero.
     */
    class halfSpherePartition : public spherePartition {
    public:
        explicit halfSpherePartition(int numSections);
    };
}