//
// Created by maojrs on 4/6/22.
//

#pragma once
#include "potentials.hpp"

namespace msmrd {
    /**
     * Class to combine several pair potentials into one. Potentials need to be added with the addPotential
     * function.
     */
    class combinedPairPotential : public pairPotential{
    public:
        std::vector<pairPotential*> potentials;
        /*
         * @param pairPotentials, tuple of pointers to pairPotentials. The evaluate and forceTorque
         * functions call all the potentials in the tuple to evaluate the functions.
         */

        using pairPotential::pairPotential;

        void addPotential(pairPotential *pairPotentialPtr);

        double evaluate(particle &part1, particle &part2) override;

        std::array<vec3<double>, 4> forceTorque(particle &part1, particle &part2) override;
    };

    /**
    * Class to combine several external potentials into one. Potentials need to be added with the
     * addPotential function.
    */
    class combinedExternalPotential : public externalPotential{
    public:
        std::vector<externalPotential*> potentials;
        /*
         * @param pairPotentials, tuple of pointers to pairPotentials. The evaluate and forceTorque
         * functions call all the potentials in the tuple to evaluate the functions.
         */

        using externalPotential::externalPotential;

        void addPotential(externalPotential *externalPotentialPtr);

        double evaluate(particle &part) override;

        std::array<vec3<double>, 2> forceTorque(particle &part) override;
    };


}