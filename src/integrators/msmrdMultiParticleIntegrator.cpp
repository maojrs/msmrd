//
// Created by maojrs on 10/7/19.
//

#include "integrators/msmrdMultiParticleIntegrator.hpp"

namespace msmrd {
    /**
     * Implementation of multiparticle MSM/RD integrator. Based mainly on msmrdintegratorDiscrete.
     * Uses same constructor as msmrdIntegratorDiscrete, with the followin parameters
     * @param dt time step
     * @param seed random generator seed (Note seed = -1 corresponds to random device)
     * @param rotation boolean to indicate if rotational degrees of freedom should be integrated
     */

    template<>
    void msmrdMultiParticleIntegrator<ctmsm>::transition2BoundState(std::vector<particleMS> &parts, int iIndex,
                                                             int jIndex, int endState) {

    }
}