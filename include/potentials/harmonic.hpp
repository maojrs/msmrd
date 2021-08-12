//
// Created by maojrs on 8/12/21.
//

#pragma once
#include "vec3.hpp"
#include "potentials/potentials.hpp"

namespace msmrd {
    /*
     * Harmonic external potential in 3D corresponding to - x*x. One can provide a coefficient in each
     * dimension for non-isotropic harmonic potential and a location for the minima.
     */
    class harmonic : public externalPotential {
    public:
        double scalefactor = 1;
        long seed = 0;
        std::vector<int> particleTypes;
        vec3<double> minima;
        vec3<double> kconstants;

        /**
         * @param scalefactor adjusts strength of the potential muliplying by a constant
         * @param seed variable for random number generation (Note seed = -1 corresponds to random device)
         * @param particleTypes if assigned values, the external potential only acts on the specified
         * particle types.
         * @param minima coordinates of the minima
         * @param kconstants 3d vectors of k constant in each dimension
         */
        harmonic(std::vector<double> minima, std::vector<double> kconstants, double scalefactor);

        harmonic(std::vector<double> minima, std::vector<double> kconstants,
                 std::vector<int> particleTypes, double scalefactor);

        double evaluate(particle &part) override;

        std::array<vec3<double>, 2> forceTorque(particle &part) override;
    };

}
