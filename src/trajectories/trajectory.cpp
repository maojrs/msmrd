//
// Created by dibakma on 18.09.18.
//

#include <iostream>
#include "trajectories/trajectory.hpp"
#include "particle.hpp"


// Needed to write to HDF5 files.
using namespace H5;

namespace msmrd {

    /**
     * Implementation of abstract parent trajectory class
     * @param Nparticles number of particles in the simulation that need to be saved in the trajectory
     * @param bufferSize buffer size (per particle) for the trajectoryData array. If not using H5, it will be
     * used as first estimate to initialize the trajectoryData array (recommended == timeIterations/stride). If
     * using H5, it will determine the size stored in memory before flushing data into file and emptyinf buffer.
     */
    trajectory::trajectory(unsigned long Nparticles, int bufferSize): Nparticles(Nparticles), bufferSize(bufferSize){};

    // Empties trajectories data buffers
    void trajectory::emptyBuffer() {
        trajectoryData.clear();
        discreteTrajectoryData.clear();
    }

    // Incorporates custom boundary into integrator
    void trajectory::setBoundary(boundary *bndry) {
        boundaryActive = true;
        domainBoundary = bndry;
    }

    // Calculates relative position using boundary information loaded into trajectory class (position2 - position1)
    vec3<double> trajectory::calculateRelativePosition(vec3<double> position1, vec3<double> position2) {
        if (boundaryActive and domainBoundary->getBoundaryType() == "periodic") {
            auto boxsize = domainBoundary->boxsize;
            return msmrdtools::distancePeriodicBox(position1, position2, boxsize);
        } else {
            return position2 - position1;
        }
    }
}