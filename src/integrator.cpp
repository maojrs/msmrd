//
// Created by dibakma on 27.06.18.
//
#include <array>
#include <utility>
#include "integrator.hpp"
#include "particle.hpp"


odLangevin::odLangevin(double dt, long seed) : integrator(dt,seed) {};

void odLangevin::integrate(std::vector<std::shared_ptr<particle>> parts) {
    vec3<double> dr;
    double coeff;
    for (int i=0; i<parts.size(); i++) {
        coeff = std::sqrt(2*dt*parts[i]->D);
        dr[0] = coeff*randg.normal(0,1);
        dr[1] = coeff*randg.normal(0,1);
        dr[2] = coeff*randg.normal(0,1);
        parts[i]->position += dr;
    }
};

double odLangevin::test(){
    return randg.normal(0,1);
}

