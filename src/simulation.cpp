//
// Created by dibakma on 22.06.18.
//

//#include <random>
#include "simulation.hpp"
#include "trajectory.hpp"


namespace msmrd {
    void simulation::run(const int Nsteps, trajectory& traj, int stride) {
        for (int step=0; step < Nsteps; step++) {
            integ.integrate(particleList);
            if (step % stride == 0) {
                traj.sample(step, particleList);
            }
        }
    }
}