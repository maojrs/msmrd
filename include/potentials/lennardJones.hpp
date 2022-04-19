//
// Created by maojrs on 4/19/21.
//

#pragma once
#include "potentials.hpp"

namespace msmrd {
    /*
     * Declaration of the Lennard Jones pair potential. We also include the WCA potential here, which is simply a
     * Lennard-Jones pontetical truncated at its minima, so it only has the repulsive part.
     */
    class lennardJones : public pairPotential {
    public:
        double epsilon;
        double sigma;
        double cutOff;
        bool forceCap = false;
        double forceCapValue = 0.0;
        double potentialCutOff = 0.0;
        double rcritical = 0;
        double baseEnergy = 0.0;
        std::vector<int> excludeParticleTypesPairs;

        /**
         * @param epsilon depth of the potential well
         * @param sigma distance at which the potential is zero (a.k.a 'size of the particle')
         * @param cutOff value after which potential is zero, defaults to 3 * sigma if not specified
         * @param forceCap establish a maximum possible value for the force regardless of potential
         * @param forceCapValue maximum possible value of force (if forceCap = true)
         * @param potentialCutOff value at which LJ potential is cut off and substituted by a line with constant
         * slope(force) (only if forceCap = true)
         * @param rcritical r value at which the LJ potential is cutoff (only if forceCap = true)
         * @param baseEnergy translation of potential energy in y-axis. It defaults to zero, but it can
         * be modified by child classes, like the WCA.
         * @param excludeParticleTypesPairs if set, the pair particles of this type (both the same), will
         * be excluded from the interaction.
         */
        lennardJones(double epsilon, double sigma);

        lennardJones(double epsilon, double sigma, double cutOff);

        lennardJones(double epsilon, double sigma, std::vector<int> exclParticleTypesPairs);

        double evaluate(particle &part1, particle &part2) override;

        std::array<vec3<double>, 4>
        forceTorque(particle &part1, particle &part2) override;

        void setForceCapValue(double forceCapVal);

        void setPotentialCutOff(double potCutOff);

        double getPotentialCutOff() { return potentialCutOff; }

        double getForceCapValue() { return forceCapValue; }
    };

    class WCA : public lennardJones {
    public:
        WCA(double epsilon, double sigma);

        WCA(double epsilon, double sigma, std::vector<int> exclParticleTypesPairs);

    };

}
