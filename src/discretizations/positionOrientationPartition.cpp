//
// Created by maojrs on 4/1/19.
//

#include "discretizations/positionOrientationPartition.hpp"

namespace msmrd {

    /*
     * Constructor creates a surface partition of the unit 3D sphere to discretize the direction of the
     * relative position, plus a volumetric partition of the unit sphere in 3D that discretizes the relative
     * orientation (unit quaternion). Together they discretize the relative position and relative
     * orientation phase space. Note the magnitude of the relative position is only discretized by the relative
     * distance cut off value.
     */
    positionOrientationPartition::positionOrientationPartition(double relativeDistanceCutOff,
                                                               int numSphericalSectionsPos,
                                                               int numRadialSectionsQuat,
                                                               int numSphericalSectionsQuat) :
            numSphericalSectionsPos(numSphericalSectionsPos),
            relativeDistanceCutOff(relativeDistanceCutOff),
            numRadialSectionsQuat(numRadialSectionsQuat),
            numSphericalSectionsQuat(numSphericalSectionsQuat) {

        sphericalPartition = std::make_unique<spherePartition>(numSphericalSectionsPos);
        quatPartition = std::make_unique<quaternionPartition>(numRadialSectionsQuat, numSphericalSectionsQuat);
        numTotalQuatSections = quatPartition->numTotalSections;
        numTotalSections = numSphericalSectionsPos * quatPartition->numTotalSections;

//        /* Obtains a good naming scheme for the states, while keeping it an integer (convenient for pyemma)
//         * numDigits is the number of digits used by numTotalSections in quaternion partition.
//         * sectionNumber = sectionNumberingMultiplier*secNumRelativePos + secNumRelativeQuat.
//         * This way there is always at least zero between secNumRelativePos  and secNumRelativeQuat*/
//        int numDigits = static_cast<int>(log10(quatPartition->numTotalSections)) + 1;
//        sectionNumberingMultiplier = std::pow(10, numDigits + 1);
    };


    /* Gets the section number in this partition given a relative position and relative orientation.
     * Note in order to get consistent results relativePosition should be rotated by particle1.quat.conjugated(),
     * which is passed in as quaternionReference. This is equivalent to fixing the spherical surface partition
     * to particle1, so they rotate together. Alternatively, the quaternionReference can be set to identity if the values are already in
     * the correct frame of reference. */
    int positionOrientationPartition::getSectionNumber(vec3<double> relativePosition,
                                                       quaternion<double> relativeOrientation,
                                                       quaternion<double> quaternionReference) {
        if (relativePosition.norm() > relativeDistanceCutOff) {
            return 0;
        }
        // Rotate relative position back to original frame of reference of particle 1
        vec3<double> fixedRelativePosition = msmrdtools::rotateVec(relativePosition, quaternionReference.conj());
        // Calculate section numbers given by sphere partition and quaternion partition.
        int secNumRelativePos = sphericalPartition->getSectionNumber(fixedRelativePosition);
        int secNumRelativeQuat = quatPartition->getSectionNumber(relativeOrientation);
        // Generate section number for positionOrientationPartition from these previous section numbers
        int sectionNumber = (secNumRelativePos - 1) * quatPartition->numTotalSections + secNumRelativeQuat;
        return sectionNumber;
    };


    /* Adds an offset on the thetas of the spherical partition and of the quaternions partition. Must
     * be a positive value and preferably smaller than all the possible dthetas */
    void positionOrientationPartition::setThetasOffset(double offset) {
        sphericalPartition->setThetasOffset(offset);
        quatPartition->setThetasOffset(offset);
    };



    /* Gets other two section numbers corresponding to the positionOrientationPartition section number.
     * The function returns a tuple (secNumRelativePos, secNumRelativeQuat): the section number corresponding
     * to the surface spherical partition, and the section number corresponding to the quaternion partition,
     * respectively. The exact intervals can then be extracted by calling
     * spherePartition->getAngles(secNumRelativePos) and quatPartition->getSectionIntervals(secNumRelativeQuat),
     * as done in getSectionIntervals. */
    std::tuple<int, int> positionOrientationPartition::getSectionNumbers(int secNumber){
        int secNumRelativePos = static_cast<int>(std::ceil(1.0*secNumber/numTotalQuatSections));
        int secNumRelativeQuat = secNumber % numTotalQuatSections + 1;
        return std::make_tuple(secNumRelativePos, secNumRelativeQuat);
    };


    /* Gets the exact intervals for all the involved partitions corresponding to a given section number. It uses
     * the function getSectionNumbers to obtain the corresponding section numbers in the spherePartition and in
     * the quaternionPartition and calls spherePartition->getAngles(secNumRelativePos) and
     * quatPartition->getSectionIntervals(secNumRelativeQuat). The first two returned values in the tuple
     * correspond to the angles intervals in the spherePartition and the next three correspond to the section
     * intervals of the quaternion partition. */
    std::tuple<std::array<double, 2>, std::array<double, 2>,
            std::array<double, 2>, std::array<double, 2>, std::array<double, 2>>
    positionOrientationPartition::getSectionIntervals(int secNumber){
        auto sectionNumbers = getSectionNumbers(secNumber);
        int secNumRelativePos = std::get<0>(sectionNumbers);
        int secNumRelativeQuat = std::get<1>(sectionNumbers);
        auto angles = sphericalPartition->getAngles(secNumRelativePos);
        auto intervals = quatPartition->getSectionIntervals(secNumRelativeQuat);
        return std::tuple_cat(angles, std::move(intervals));
    };

    // Gets section number version PyBind
    int positionOrientationPartition::getSectionNumberPyBind(std::vector<double> relpos,
                                                             std::vector<double> relquat, std::vector<double> qref ){
        vec3<double> relposition{relpos[0], relpos[1], relpos[2]};
        quaternion<double> relquaternion{relquat[0], relquat[1], relquat[2], relquat[3]};
        quaternion<double> quaternionref{qref[0], qref[1], qref[2], qref[3]};
        return getSectionNumber(relposition, relquaternion, quaternionref);
    }

    // Returns total number of sections, num radial sections and number of spherical sections (for pybind)
    std::vector<int> positionOrientationPartition::getNumSections(){
        return std::vector<int>{numTotalSections, numSphericalSectionsPos,
                                numRadialSectionsQuat, numSphericalSectionsQuat};
    }


}
