//
// Created by maojrs on 3/28/19.
//
#pragma once
#include "spherePartition.hpp"

namespace msmrd {
    /*
     * This class creates an equal area partition on the surface of a sphere. It is a c++ copy and extension
     * of the python code in module msmrd2.tools.spherePartition.
     */
    class halfSpherePartition : public spherePartition {
    private:
        double scaling = 2; // Scales to half a sphere, other scalings must be carefully handled
    public:
        explicit halfSpherePartition(int numSections);

        // Given a section number give, phi and theta angles that define the section.
        std::tuple<std::vector<double>, std::vector<double>> getAnglesHalf(int secNumber);


    };
}