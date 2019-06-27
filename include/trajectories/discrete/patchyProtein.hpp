//
// Created by maojrs on 4/8/19.
//

#pragma once
#include <memory>
#include "trajectoryPositionOrientation.hpp"
#include "discretizations/positionOrientationPartition.hpp"
#include "tools.hpp"

namespace msmrd {
    /**
     * Trajectory class for patchy protein trajectories. Patchy protein here refers to two
     * patchy particles each with a different number of patches. The main difference with the general class
     * trajectoryPositionOrientation is that this class can sample the discrete trajectory that is
     * specific for the patchy protein application. In short words,  it chooses how to discretize the full
     * trajectory of two particles into a discretized trajectory to be analyzed and extracted into
     * a Markov state model.
     *
     * Patchy protein will have one unbound state (0), two (or more) bound states (1,2) and
     * Note in the current indexing implementation only up to 9 bound states are supported.
     */
    class patchyProtein : public trajectoryPositionOrientation {
    private:
        std::vector<quaternion<double>> rotMetastableStates;
        std::unique_ptr<positionOrientationPartition> positionOrientationPart;
        double tolerance = 0.15;
        int prevsample = 0;
        /*
         * @param rotMetastableStates[X] correspond to the list of relative rotations that correspond to a
         * metastable state/region X. Each rotation is represented by a quaternion.
         * @param positionOrientationPart pointer to partition class to discretize relative position and orientation.
         * @param tolerance is the maximum distance away from metastable state to still be considered metastable.
         * @param prevsample keeps calue of previous sample when sampling discrete trajectory, useful for
         * coreMSM approach
         */
    public:


        patchyProtein(unsigned long Nparticles, int bufferSize,
                      double relativeDistanceCutOff, int numSphericalSectionsPos,
                      int numRadialSectionsQuat, int numSphericalSectionsQuat);

        void sampleDiscreteTrajectory(double time, std::vector<particle> &particleList) override;

        void setMetastableRegions();

    };


}
