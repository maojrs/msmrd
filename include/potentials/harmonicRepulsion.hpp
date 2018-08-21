//
// Created by maojrs on 8/16/18.
//

#pragma once
#include "vec3.hpp"
#include "potentials/potential.hpp"

/*
 * Harmonic repulsion pair potential definition
 */
class harmonicRepulsion: public pairPotential{
public:
    double k;
    double range;
    /**
     * @param k repulsion strength
     * @param range interaction radius
     */
    harmonicRepulsion(double k, double range);

    double evaluate(vec3<double> pos1, vec3<double> pos2) override;
    double evaluatePyBind(std::vector<double> pos1, std::vector<double> pos2);
    vec3<double> force(vec3<double> pos1, vec3<double> pos2) override;
    std::vector<double> forcePyBind(std::vector<double> pos1, std::vector<double> pos2);
};