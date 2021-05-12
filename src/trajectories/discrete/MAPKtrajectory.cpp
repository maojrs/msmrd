//
// Created by maojrs on 5/3/21.
//

#include "trajectories/discrete/MAPKtrajectory.hpp"

namespace msmrd {

    /*
     * Constructors and definition of specialized MAPK trajectory
     * IMPORTANT: Need to define pointer to positionOrientvectorPartition and
     * call setBoundStates in constructors of childs of discreteTrajectory.
     *
     * Note it uses several different varaiables and fucntions than the parent class,
     * discreteTrajectory. Be careful when calling functions from the parent class as
     * functionality might break.
     */

    MAPKtrajectory::MAPKtrajectory(unsigned long Nparticles, int bufferSize) :
            discreteTrajectory(Nparticles, bufferSize) {

        setRadialBounds(1.25, 2.25);
        setTolerances(0.12, 0.12);
        anglePatches = M_PI/2;
        int numSphericalSectionsPos = 6;
        int numSphericalSectionsOrientvec = 6;
        positionOrientvectorPart = std::make_unique<positionOrientvectorPartition>(rUpperBound,
                                                                                   numSphericalSectionsPos, numSphericalSectionsOrientvec);

        // Set bound states defined for the patchy protein example
        setBoundStates();
    };

    MAPKtrajectory::MAPKtrajectory(unsigned long Nparticles, int bufferSize, double anglePatches) :
            discreteTrajectory(Nparticles, bufferSize), anglePatches(anglePatches) {

        setRadialBounds(1.25, 2.25);
        setTolerances(0.12, 0.12);

        int numSphericalSectionsPos = 6;
        int numSphericalSectionsOrientvec = 6;
        positionOrientvectorPart = std::make_unique<positionOrientvectorPartition>(rUpperBound,
                numSphericalSectionsPos, numSphericalSectionsOrientvec);

        // Set bound states defined for the patchy protein example
        setBoundStates();
    };

    MAPKtrajectory::MAPKtrajectory(unsigned long Nparticles, int bufferSize, double anglePatches,
            double rLowerBound, double rUpperBound) :
            discreteTrajectory(Nparticles, bufferSize, rLowerBound, rUpperBound),
            anglePatches(anglePatches) {

        setTolerances(0.12, 0.12);

        int numSphericalSectionsPos = 6;
        int numSphericalSectionsOrientvec = 6;
        positionOrientvectorPart = std::make_unique<positionOrientvectorPartition>(rUpperBound,
                numSphericalSectionsPos, numSphericalSectionsOrientvec);

        // Set bound states defined for the patchy protein example
        setBoundStates();
    };

    MAPKtrajectory::MAPKtrajectory(unsigned long Nparticles, int bufferSize, double anglePatches,
            int numSphericalSectionsPos, int numSphericalSectionsOrientvec,
            double rLowerBound, double rUpperBound) :
            discreteTrajectory(Nparticles, bufferSize, rLowerBound, rUpperBound),
            anglePatches(anglePatches) {

        setTolerances(0.12, 0.12);

        positionOrientvectorPart = std::make_unique<positionOrientvectorPartition>(rUpperBound,
                numSphericalSectionsPos, numSphericalSectionsOrientvec);

        // Set bound states defined for the patchy protein example
        setBoundStates();
    };



    /* Sets bound states (metastable regions) of the MAPK implementation. The centers of the
    * metastable regions are given by a tuple of relative position, relative orientvector and type.
     * The size of the regions are determined by tolerancePosition and toleranceOrientation.
     * The relative orientvector is simply the orientvector of particle 2, seen from the frame
     * of reference of particle 1. Note here particle 1 will always be the main MAPK particle
     * and particle 2 will always be either the kinase or the phosphatase. */
    void MAPKtrajectory::setBoundStates() {
        /* Define relative position vectors from particle 1 at the origin. These two patches
         * point in the same direction as the two patches in the MAPK. */
        vec3<double> relPos1 = {std::cos(anglePatches / 2.0), std::sin(anglePatches / 2.0), 0};
        vec3<double> relPos2 = {std::cos(anglePatches / 2.0), std::sin(-anglePatches / 2.0), 0};
        /* Relative orientation vectors (assuming particle 1 --MAPK-- fixed) of particle 2 that
         * yield the bound states. */
        std::array<vec3<double>, 2> orientVecs;
        orientVecs[0] = -1 * relPos1; //bound in patch1
        orientVecs[1] = -1 * relPos2; //bound in patch2

        /* Fill bound states with corresponding combinations of relative position vectors and quaternion orientations.
         * Note we define relativeOrientation as q2 * q1.conj(), so it matches the relative orientation from
         * particle 1 as used in trajectories/discrete/discreteTrajectory.hpp */
        boundStates[0] = std::make_tuple(relPos1, orientVecs[0], 1); // kinase (type 1) in patch1
        boundStates[1] = std::make_tuple(relPos2, orientVecs[1], 1); // kinase (type 1) in patch2
        boundStates[2] = std::make_tuple(relPos1, orientVecs[0], 2); // phosph (type 2) in patch1
        boundStates[3] = std::make_tuple(relPos2, orientVecs[1], 2); // phosph (type 2) in patch2
    }

    /* Main function to sample the discrete state of two particles. It returns the corresponding
     * bound state, transition state or unbound state (0). In the bound region (r< rLowerBound), it can
     * also return -1 when not in any bound state. In this case, one would normally apply the coreMSM
     * approach and choose the previous value. However, this is done directly on sampleDiscreteTrajectory or
     * in discretizeTrajectoryH5 and discretizeTrajectory if discretizing directly a python array. */
    int MAPKtrajectory::sampleDiscreteState(particle partA, particle partB) {
        // Enforce that particle one is always the MAPK and particle 2 always the kinase or phosphastase
        particle* part1;
        particle* part2;
        if (partA.type == 0 and partB.type > 0) {
            part1 = &partA;
            part2 = &partB;
        }
        else if (partB.type == 0 and partA.type > 0) {
            part1 = &partB;
            part2 = &partA;
        }
        else{
            return 0;
        }

        // Initialize sample with value zero (unbound state)
        int discreteState = 0;

        /* Calculate relative position taking into account periodic boundary measured
         * from i to j (gets you from i to j). */
        auto relativePosition = calculateRelativePosition(part1->position, part2->position);
        /* Rotate relative position and orientVector to match the reference orientation of
         * particle 1. (VERY IMPORTANT) */
        relativePosition = msmrdtools::rotateVec(relativePosition, part1->orientation.conj());
        auto relativeOrientvector = msmrdtools::rotateVecOffAxis(part1->orientvector,
                part1->orientation.conj(), relativePosition);

        // Extract current state, save into sample and return sample
        if (relativePosition.norm() < rLowerBound) {
            // Returns discrete state or -1 if it is not in any bound state
            discreteState = getBoundState(relativePosition, relativeOrientvector, part2->type);
        }
            // Returns a transitions state if it is in the transition region
        else if (relativePosition.norm() < positionOrientvectorPart->relativeDistanceCutOff) {
            // Get corresponding section numbers from spherical partition to classify its state
            auto secNum = positionOrientvectorPart->getSectionNumber(relativePosition,relativeOrientvector);
            discreteState  = maxNumberBoundStates + secNum;
        }
        // If none of the statements before modified discreteState, it returns the unbound state (0)
        return discreteState;
    }

    int MAPKtrajectory::getBoundState(vec3<double> relativePosition, vec3<double> orientVector,
            int ligandType) {
        // Check if it matches a bound states, if so return the corresponding state, otherwise return -1.
        for (int i = 0; i < boundStates.size(); i++) {
            auto relPosCenter = std::get<0>(boundStates[i]);
            auto relOrientvec = std::get<1>(boundStates[i]);
            auto ligType = std::get<2>(boundStates[i]);
            if ( (relPosCenter - relativePosition).norm() <= tolerancePosition) {
                if  ((relOrientvec - orientVector).norm() <= toleranceOrientation) {
                    if (ligType == ligandType) {
                        return i + 1;
                    }
                }
            }
        }
        return -1;
    };

    vec3<double> MAPKtrajectory::getRelativePosition(int boundStateIndex){
        return std::get<0>(boundStates[boundStateIndex]);
    };

    quaternion<double> MAPKtrajectory::getRelativeOrientvector(int boundStateIndex){
        return std::get<1>(boundStates[boundStateIndex]);
    };

    void MAPKtrajectory::setAnglePatches(double newAnglePatches) {
        anglePatches = newAnglePatches;
        setBoundStates();
    };


}
