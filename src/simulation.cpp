//
// Created by dibakma on 22.06.18.
//

//#include <random>
#include "simulation.hpp"

template<typename rng>
int something() {
    std::mt19937 generator;
    generator.seed(10);
    std::normal_distribution<double> distribution(0.0, 1.0);
    distribution(generator);
}

//void simulation::run(const double timestep, const int Nsteps) {
//    for (int step; step<Nsteps; step++) {
//        for (int i; i<Nparticles; i++) {
//            particles.at(i).position += 1;
//        }
//    }
//}