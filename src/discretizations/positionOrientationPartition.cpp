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
        numTotalSections = numSphericalSectionsPos * quatPartition->numTotalSections;

//        /* Obtains a good naming scheme for the states, while keeping it an integer (convenient for pyemma)
//         * numDigits is the number of digits used by numTotalSections in quaternion partition.
//         * sectionNumber = sectionNumberingMultiplier*secNumRelativePos + secNumRelativeQuat.
//         * This way there is always at least zero between secNumRelativePos  and secNumRelativeQuat*/
//        int numDigits = static_cast<int>(log10(quatPartition->numTotalSections)) + 1;
//        sectionNumberingMultiplier = std::pow(10, numDigits + 1);
    };

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


    std::tuple<int, int> positionOrientationPartition::getSectionNumbers(int secNumber){
        int secNumRelativePos = static_cast<int>(std::ceil(1.0*secNumber/numSphericalSectionsQuat));
        int secNumRelativeQuat = secNumber % numSphericalSectionsQuat + 1;
        return std::make_tuple(secNumRelativePos, secNumRelativeQuat);
    };

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


}