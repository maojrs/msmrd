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
        int numSphericalSectionsOrientvec = 6;
        positionOrientvectorPart = std::make_unique<positionOrientvectorPartition>(rUpperBound,
                numSphericalSectionsPos, numSphericalSectionsOrientvec);

        // Set bound states defined for the patchy protein example
        setBoundStates();
    };

    MAPKtrajectory::MAPKtrajectory(unsigned long Nparticles, int bufferSize,
            double rLowerBound, double rUpperBound) :
            discreteTrajectory(Nparticles, bufferSize, rLowerBound, rUpperBound) {

        setTolerances(0.12, 0.12*2*M_PI);

        int numSphericalSectionsPos = 6;
        int numSphericalSectionsOrientvec = 6;
        positionOrientvectorPart = std::make_unique<positionOrientvectorPartition>(rUpperBound,
                numSphericalSectionsPos, numSphericalSectionsOrientvec);

        // Set bound states defined for the patchy protein example
        setBoundStates();
    };

    void MAPKtrajectory::setBoundStates() {

    }

    int MAPKtrajectory::sampleDiscreteState(particle part1, particle part2) {
        return 0;
    }

    int MAPKtrajectory::getBoundState(vec3<double> relativePosition, vec3<double> orientVector) {
        return 0;
    };

    vec3<double> MAPKtrajectory::getRelativePosition(int boundStateIndex){
        return vec3<double>(0,0,0);
    };

    quaternion<double> MAPKtrajectory::getRelativeOrientation(int boundStateIndex){
        return quaternion<double>();
    };


}
