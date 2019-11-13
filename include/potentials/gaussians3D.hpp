//
// Created by maojrs on 8/16/18.
//

#pragma once
#include "vec3.hpp"
#include "potentials/potentials.hpp"

namespace msmrd {
    /*
     * 3D external potential composed of nminima Gaussians placed randomly inside sphere of radius maxrad.
     * The empty template <> indicates orientation is not taken into account by this potential.
     */
    class gaussians3D : public externalPotential {
    private:
        randomgen randg;
    public:
        int nminima;
        double maxrad;
        double scalefactor;
        long seed;
        std::vector<vec3<double>> minimas;
        std::vector<vec3<double>> stddevs;

        /**
         * @param nminima number of minimas/gaussians in the potential
         * @param maxrad maximum radius where all centers of Gaussians will be contained
         * @param scalefactor adjusts strength of the potential muliplying by a constant
         * @param seed variable for random number generation (Note seed = -1 corresponds to random device)
         * @param minimas stores array of coordinates of the minimas/centers of the Gaussians
         * @param stddevs stores array of 3d vectors with the standard deviation of the Gaussian in each dimension
         */
        gaussians3D(unsigned long nminima, double maxrad, double scalefactor, long seed);

        double evaluate(particle &part) override;

        std::array<vec3<double>, 2> forceTorque(particle &part) override;
    };

}