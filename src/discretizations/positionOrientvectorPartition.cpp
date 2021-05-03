//
// Created by maojrs on 5/3/21.
//

#include "discretizations/positionOrientvectorPartition.hpp"

namespace msmrd {

    /*
     * Constructor creates a surface partition of the unit 3D sphere to discretize the direction of the
     * relative position, plus another surface partition of the unit sphere that discretizes the relative
     * orientation (orientvector, 2 degrees of freedom only). Together they discretize the relative
     * position and relative orientation phase space. Note the magnitude of the relative position is
     * only discretized by the relative distance cut off value.
     */
    positionOrientvectorPartition::positionOrientvectorPartition(double relativeDistanceCutOff,
            int numSphericalSectionsPos,int numSphericalSectionsOrientvec) :
            relativeDistanceCutOff(relativeDistanceCutOff),
            numSphericalSectionsPos(numSphericalSectionsPos),
            numSphericalSectionsOrientvec(numSphericalSectionsOrientvec){

        sphericalPartition = std::make_unique<spherePartition>(numSphericalSectionsPos);
        sphericalPartitionOrientvec = std::make_unique<spherePartition>(numSphericalSectionsOrientvec);
        numTotalSections = numSphericalSectionsPos * numSphericalSectionsOrientvec;
    };


    /* Gets the section number in this partition given a relative position and relative orientvector.
     * Note in order to get consistent results, relativePosition and relativeOrientvector should be
     * rotated by particle1.quat.conjugated() (relativeOrientvector off axis rotation). This is
     * equivalent to fixing the spherical surface partition to particle1, so they rotate together. The
     * function assumes the relativePosition and orientvector are already in the frame of reference
     * fixed to the main particle.*/
    int positionOrientvectorPartition::getSectionNumber(vec3<double> relativePosition,
                                                       vec3<double> orientvector) {
        if (relativePosition.norm() > relativeDistanceCutOff) {
            return 0;
        }
        // Calculate section numbers given by sphere partition and quaternion partition.
        int secNumRelativePos = sphericalPartition->getSectionNumber(relativePosition);
        int secNumRelativeQuat = sphericalPartitionOrientvec->getSectionNumber(orientvector);
        // Generate section number for positionOrientationPartition from these previous section numbers
        int sectionNumber = (secNumRelativePos - 1) * numSphericalSectionsOrientvec + secNumRelativeQuat;
        return sectionNumber;
    };


    /* Adds an offset on the thetas of both spherical partitions. Must
     * be a positive value and preferably smaller than all the possible dthetas */
    void positionOrientvectorPartition::setThetasOffset(double offset) {
        sphericalPartition->setThetasOffset(offset);
        sphericalPartitionOrientvec->setThetasOffset(offset);
    };



    /* Gets other two section numbers corresponding to the positionOrientvectorPartition section number.
     * The function returns a tuple (secNumRelativePos, secNumRelativeOrientvec): the section number
     * corresponding to the surface spherical partition of the position, and the section number
     * corresponding to the surface spherical partition for the orientvector, respectively. The exact
     * intervals can then be extracted by calling spherePartition->getAngles() as done in
     * getSectionIntervals. */
    std::tuple<int, int> positionOrientvectorPartition::getSectionNumbers(int secNumber){
        int secNumRelativePos = static_cast<int>(std::ceil(1.0*secNumber/numSphericalSectionsOrientvec));
        int secNumRelativeQuat = secNumber % numSphericalSectionsOrientvec + 1;
        return std::make_tuple(secNumRelativePos, secNumRelativeQuat);
    };


    /* Gets the exact intervals for all the involved partitions corresponding to a given section number. It uses
     * the function getSectionNumbers to obtain the corresponding section numbers in the spherePartitions
     * and calls spherePartition->getAngles. The first two returned values in the tuple
     * correspond to the angles intervals in the spherePartition and the other two correspond to the section
     * angles interval for the orientvector. */
    std::tuple<std::array<double, 2>, std::array<double, 2>,
            std::array<double, 2>, std::array<double, 2>>
    positionOrientvectorPartition::getSectionIntervals(int secNumber){
        auto sectionNumbers = getSectionNumbers(secNumber);
        int secNumRelativePos = std::get<0>(sectionNumbers);
        int secNumRelativeOrientvec = std::get<1>(sectionNumbers);
        auto angles = sphericalPartition->getAngles(secNumRelativePos);
        auto anglesOrientvec = sphericalPartitionOrientvec->getAngles(secNumRelativeOrientvec);
        return std::tuple_cat(angles, std::move(anglesOrientvec));
    };

    // Gets section number version PyBind
    int positionOrientvectorPartition::getSectionNumberPyBind(vec3<double> relpos,
                                                              vec3<double> orientvec){
        vec3<double> relposition{relpos[0], relpos[1], relpos[2]};
        vec3<double> relorientvec{orientvec[0], orientvec[1], orientvec[2]};
        return getSectionNumber(relposition, relorientvec);
    }

    // Returns total number of sections, num radial sections and number of spherical sections (for pybind)
    std::vector<int> positionOrientvectorPartition::getNumSections(){
        return std::vector<int>{numTotalSections, numSphericalSectionsPos, numSphericalSectionsOrientvec};
    }


}
