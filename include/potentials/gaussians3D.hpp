//
// Created by maojrs on 8/16/18.
//

#pragma once
#include "vec3.hpp"
#include "potentials/potential.hpp"

/*
 * 3D external potential composed of nminima Gaussians placed randomly inside sphere of radius maxrad
 */
class gaussians3D: public potential {
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
    gaussians3D(int nminima, double maxrad, double scalefactor, long seed);

    // Functions with PyBind suffix required for Pybind (compatible with std::vector)
    double evaluate(vec3<double> pos) override;
    double evaluatePyBind(std::vector<double> pos);
    vec3<double> force(vec3<double> pos) override;
    std::vector<double> forcePyBind(std::vector<double> pos);
};