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

        // Wrappers for parent functions:

        /* Given a section number give, phi and theta angles that define the section.
         * Note it may hide the parent function. */
        std::tuple<std::vector<double>, std::vector<double>> getAngles(int secNumber);

        /* Given a coordinate return corresponding section number.
         * Note it may hide the parent function. */
        int getSectionNumber(vec3<double> coordinate);


    };
}