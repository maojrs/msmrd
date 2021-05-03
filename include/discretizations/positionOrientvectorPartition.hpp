//
// Created by maojrs on 5/3/21.
//

#pragma once
#include <memory>
#include "spherePartition.hpp"
#include "tools.hpp"


namespace msmrd {
    /*
     * This class discretizes the relative position and relative orientation (given by an orientation vector, i.e.
     * only two degrees of freedom) of two rigid bodies. This can be used to create MSM/RD discretizations
     * for a generic rigid body particle (orientation given by quaternion) interacting with a rod-like particle
     * (orientation given by orientvector). To do the discretization, it uses the spherePartition twice. The
     * first one to discretize the relative position vector (mainly its direction), and the latter one to
     * discretize the relative orientation. Overall the space spanned by the relative position and orientation
     * is effectively 5 dimensional, so this provides effectively a discretization of a 5 dimensional
     * space. This is a somewhat specialized class aimed for specific applications. For generic applications
     * stick to positionOrientationPartition.
     *
     * Note the surface sphere partition to discretize the relative position should be fixed to
     * one of the particles, this should be the generic rigid body particle (orientation described
     * by quaternion). The partition rotates with this particle. As it is simpler to leave the
     * partition in place, we simply rotate the relative position back to the initial orientation
     * of particle 1.
     *
     * Note this class in general does need to know the orientation of particles in terms of quaternions, so
     * it should be implemented with rigid bodies. However, to determine the discretization it only
     * uses the orientvector, which should be updated simoultaneously to the orientation by the integrator.
     *
     * Note section numbering (secNumber) starts in 1 and not zero.
     */
    class positionOrientvectorPartition {
    public:
        double relativeDistanceCutOff;
        int numSphericalSectionsPos;
        int numSphericalSectionsOrientvec;
        int numTotalSections = 0;
        std::unique_ptr<spherePartition> sphericalPartition;
        std::unique_ptr<spherePartition> sphericalPartitionOrientvec;
        /**
         * @param reltiveDistanceCutOff determines the cutoff for the relative position, after this value is exceeded,
         * the partition is no longer effective and the getSectionNumber function returns 0 (unbound state)
         * @param numSphericalSectionsPos number of sections in the surface partition of the unit sphere to discretize
         * the direction of the relative position.
         * @param numSphericalSectionsOrientvec number of sections in the surface partition of the unit sphere to
         * discretize the direction of the relative orientation (simple unit vector).
         * @param numTotalSections total number of sections used to discretize relative position and orientation
         * @param sphericalPartition pointer to spherical partition that discretizes the direction
         * of the relative position vector in the surface of a unit 3D sphere.
         * @param sphericalPartitionOrientvec pointer to sphere partition that discretizes the relative
         * orientation (orientvector) in the area of a unit sphere.
         */

        positionOrientvectorPartition(double relativeDistanceCutOff, int numSphericalSectionsPos,
                                     int numSphericalSectionsOrientvec);

        void setThetasOffset(double offset);

        int getSectionNumber(vec3<double> relativePosition, vec3<double> orientvector);

        std::tuple<int, int> getSectionNumbers(int secNumber);

        std::tuple<std::array<double, 2>, std::array<double, 2>,
                std::array<double, 2>, std::array<double, 2>>
        getSectionIntervals(int secNumber);


        /* Other not so important functions (mostly for PyBindings)*/

        std::vector<int> getNumSections();

        int getSectionNumberPyBind(vec3<double> relpos, vec3<double> orientvec);

    };

}
