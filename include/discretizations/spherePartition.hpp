//
// Created by maojrs on 1/31/19.
//

#pragma once
#include <vector>
#include <tuple>
#include "vec3.hpp"
#include "quaternion.hpp"


namespace msmrd {
   /*
    * This class creates an equal area partition on the surface of a sphere. It is a c++ copy and extension
    * of the python code in module msmrd2.tools.spherePartition.
    *
    * Note section numbering (secNumber) starts in 1 and not zero.
    */
    class spherePartition {
    protected:

        void partitionSphere();

        /*The following three functions should remain outside of the main partition calculation,
         * so it can be easily generalizable to other discretizations (like the half sphere) */

        double angle2CapArea(double phi);

        double capArea2Angle(double area);

        double stateArea();

    public:
        int numSections;
        double scaling = 1;
        double thetasOffset = 0.0;
        std::vector<int> regionsPerCollar;
        std::vector<double> phis;
        std::vector<std::vector<double>> thetas;
       /* @param numSections the total number of sections of equal area in the partitions
        * The partition is fully encoded into three variables within the class:
        * @param scaling can be used to scale the partition to the half sphere (scaling = 2 instead of 1). However
        * its behavior is not trivial so modify with care.
        * @param thetasOffset is the offset to the collars in theta. This could be useful to adapt the discretization
        * to a given problem. Default value is zero, but it can be changed with setOffsetThetas function.
        * @param regionsPerCollar, vector which size denotes the number of horizontal collars to
        * split the sphere, and the value of each entry the number of sections in each collar. Summing all these
        * values should result in numPartitions
        * @param phis are the locations in the polar angle that one should make the cuts to obtain the corresponding
        * collars.
        * @param thetas each entry of this vector corresponds to one collar (except for the polar caps). Each
        * entry is a another vector with the location in the azimutal angle to make the cuts to obtain the
        * sections in each collar.
        *
        */

        spherePartition() = default;

        spherePartition(int numSections);

        void setThetasOffset(double offset);

        int getSectionNumber(vec3<double> coordinate);

        std::tuple<std::array<double, 2>, std::array<double, 2>> getAngles(int secNumber);


        /* Other not so important functions (mostly for PyBindings)*/

        std::tuple<std::vector<int>, std::vector<double>, std::vector<std::vector<double>>> getPartition();

        int getNumSections(){ return numSections; }

        int getSectionNumberPyBind(std::vector<double> coord);

    };

}