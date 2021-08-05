//
// Created by maojrs on 8/16/18.
//

#pragma once
#include "vec3.hpp"
#include "potentials/potentials.hpp"

namespace msmrd {
    /*
     * 3D external potential composed of nminima Gaussians placed randomly inside sphere of radius maxrad.
     * The empty template <> indicates orientation is not taken into account by this potential. One can also
     * provide a list of minima positions and standard deviations in an alternative constructor to avoid any
     * randomness.
     */
    class gaussians3D : public externalPotential {
    private:
        randomgen randg;
    public:
        int nminima;
        double maxrad = 1;
        double scalefactor = 1;
        long seed = 0;
        std::vector<int> particleTypes;
        std::vector<vec3<double>> minimas;
        std::vector<vec3<double>> stddevs;

        /**
         * @param nminima number of minimas/gaussians in the potential
         * @param maxrad maximum radius where all centers of Gaussians will be contained
         * @param scalefactor adjusts strength of the potential muliplying by a constant
         * @param seed variable for random number generation (Note seed = -1 corresponds to random device)
         * @param particleTypes if assigned values, the external potential only acts on the specified
         * particle types.
         * @param minimas stores array of coordinates of the minimas/centers of the Gaussians
         * @param stddevs stores array of 3d vectors with the standard deviation of the Gaussian in each dimension
         */
        gaussians3D(unsigned long nminima, double maxrad, double scalefactor, long seed);

        gaussians3D(std::vector<std::vector<double>> minimaPositions,
                std::vector<std::vector<double>> standardDeviations, double scalefactor);

        gaussians3D(std::vector<std::vector<double>> minimaPositions,
                    std::vector<std::vector<double>> standardDeviations,
                    std::vector<int> particleTypes, double scalefactor);

        double evaluate(particle &part) override;

        std::array<vec3<double>, 2> forceTorque(particle &part) override;
    };

}