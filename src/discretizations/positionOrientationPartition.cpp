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
        numTotalSections = numSphericalSectionsPos*quatPartition->numTotalSections;

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

    std::tuple<int, int> positionOrientationPartition::getSection(int secNumber){
        int secNumRelativePos = static_cast<int>(std::ceil(1.0*secNumber/numSphericalSectionsQuat));
        int secNumRelativeQuat = secNumber % numSphericalSectionsQuat + 1;
        return std::make_tuple(secNumRelativePos, secNumRelativeQuat);
    };


}