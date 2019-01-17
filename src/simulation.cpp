//
// Created by dibakma on 22.06.18.
//

//#include <random>
#include "simulation.hpp"
#include <memory>
#include "trajectories/trajectory.hpp"


namespace msmrd {

    simulation::simulation(std::vector<particle> &particleList,  integrator &integ)
            : particleList(particleList), integ(integ) {};

    void simulation::run(const int Nsteps, const int buffersize, const int stride, const std::string filename) {
        int bufferCounter = 0;
        auto p1 = vec3<double> {0.0, 0.0, 0.0};
        auto p2 = vec3<double> {1.0, 0.0, 0.0};
        auto o1 = quaternion<double> {1.0, 0.0, 0.0, 0.0};
        auto o2 = quaternion<double> {1.0, 0.0, 0.0, 0.0};
        particle part1(1., 1., p1, o1);
        particle part2(1., 1., p2, o2);
        std::vector<particle> particlelist {part1, part2};


        if (integ.getParticlesBodyType() == "point"){
            traj = std::make_unique<trajectoryPosition>(particlelist.size(), buffersize);
        } else {
            traj = std::make_unique<trajectoryPositionOrientation>(particlelist.size(), buffersize);
        }

        trajectoryPositionOrientation traj2 = trajectoryPositionOrientation(particlelist.size(), buffersize);


        for (int tstep=0; tstep < Nsteps; tstep++) {
            integ.integrate(particlelist);
            if (tstep % stride == 0) {
                bufferCounter++;
                traj->sample(tstep, particlelist);
                //if (bufferCounter == traj.bufferSize) {
                //    bufferCounter = 0;

                //}
            }
        }
        traj->write2H5file<8>(filename, traj->getData());
        traj->write2file(filename, traj->getData());
    }
}