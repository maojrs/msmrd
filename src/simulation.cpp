//
// Created by dibakma on 22.06.18.
//

//#include <random>
#include "simulation.hpp"
#include "trajectories/trajectory.hpp"


namespace msmrd {

    simulation::simulation(integrator &integ, std::vector<particle> &particleList)
            : integ(integ), particleList(particleList) {};

    void simulation::run(const int Nsteps, trajectory& traj, int stride) {
        int bufferCounter = 0;
        std::string filename = "test.h5";
        for (int step=0; step < Nsteps; step++) {
            integ.integrate(particleList);
            if (step % stride == 0) {
                bufferCounter++;
                traj.sample(step, particleList);
                if (traj.bufferSize) {
                    bufferCounter = 0;
                    traj.write2H5file(filename, traj.data);
                }
            }
        }
    }
}