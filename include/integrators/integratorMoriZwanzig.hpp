//
// Created by maojrs on 6/14/21.
//

#pragma once
#include "integrators/langevin.hpp"
#include "particle.hpp"
#include "potentials/potentials.hpp"

namespace msmrd {
    /**
     * Integrator class for integratorMoriZwanzig. It is essentially the same as the langevin integrator, but
     * it stores specific values into the raux variable of the particle class for further analysis.
     *
     * It assumes the distinguished particle is the one with index0 (type 1).
     */
    class integratorMoriZwanzig : public langevin {
    protected:
        std::vector<int> distinguishedTypes{1};
        /**
         * @distinguishedTypes vector of types that correspond to distinguished particle. Must match that of trajectory.
         */

        void loadAuxiliaryValues(std::vector<particle> &parts);
    public:
        using langevin::langevin;

        void integrate(std::vector<particle> &parts) override;

        void setDistinguishedTypes(std::vector<int> newDistinguishedTypes) {
            distinguishedTypes = newDistinguishedTypes; }



    };
}