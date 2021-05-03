//
// Created by maojrs on 5/3/21.
//

#include "trajectories/discrete/MAPKtrajectory.hpp"

namespace msmrd {

    /*
     * Constructors and definition of bound states of patchy protein classes.
     * IMPORTANT: Need to define pointer to positionOrientationPartition and
     * call setBoundStates in constructors of childs of discreteTrajectory.
     */

    MAPKtrajectory::MAPKtrajectory(unsigned long Nparticles, int bufferSize) :
            discreteTrajectory(Nparticles, bufferSize) {

        setRadialBounds(1.25, 2.25);
        setTolerances(0.12, 0.12*2*M_PI);

        int numSphericalSectionsPos = 6;
        int numRadialSectionsQuat = 4;
        int numSphericalSectionsQuat = 6;
        positionOrientationPart = std::make_unique<positionOrientationPartition>(rUpperBound, numSphericalSectionsPos,
                                                                                 numRadialSectionsQuat,
                                                                                 numSphericalSectionsQuat);

        // Set bound states defined for the patchy protein example
        setBoundStates();
    };

    MAPKtrajectory::MAPKtrajectory(unsigned long Nparticles, int bufferSize,
            double rLowerBound, double rUpperBound) :
            discreteTrajectory(Nparticles, bufferSize, rLowerBound, rUpperBound) {

        setTolerances(0.12, 0.12*2*M_PI);

        int numSphericalSectionsPos = 6;
        int numRadialSectionsQuat = 4;
        int numSphericalSectionsQuat = 6;
        positionOrientationPart = std::make_unique<positionOrientationPartition>(rUpperBound, numSphericalSectionsPos,
                                                                                 numRadialSectionsQuat,
                                                                                 numSphericalSectionsQuat);

        // Set bound states defined for the patchy protein example
        setBoundStates();
    };

    void MAPKtrajectory::setBoundStates() {

    }

    int MAPKtrajectory::sampleDiscreteState(particle part1, particle part2) {
        return 0;
    }


}
