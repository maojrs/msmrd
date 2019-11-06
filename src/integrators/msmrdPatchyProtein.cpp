//
// Created by maojrs on 11/5/19.
//

#include "integrators/msmrdPatchyProtein.hpp"

namespace msmrd {
    /**
     * Implementation of MSM/RD integrator for patchyProtein integration.
     */

    int msmrdPatchyProtein::computeCurrentTransitionState(particle &part1, particle &part2) {
        int currentTransitionState = msmrdIntegrator<ctmsm>::computeCurrentTransitionState(part1, part2);
        if (part2.state == 1 and currentTransitionState > markovModel.getMaxNumberBoundStates()) {
            currentTransitionState += positionOrientationPart->numTotalSections;
        }
        return currentTransitionState;
    }

}