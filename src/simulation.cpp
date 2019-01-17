//
// Created by dibakma on 22.06.18.
//

//#include <random>
#include "simulation.hpp"
#include <memory>
#include "trajectories/trajectory.hpp"


namespace msmrd {

    simulation::simulation(std::vector<particle> &particleList,  integrator &integ)
            : particleList(particleList), integ(integ){};

    void simulation::run(const int Nsteps, const int bufferSize, const int stride, const std::string filename) {
        int bufferCounter = 0;
        auto p1 = vec3<double> {0.0, 0.0, 0.0};
        auto p2 = vec3<double> {1.0, 0.0, 0.0};
        auto o1 = quaternion<double> {1.0, 0.0, 0.0, 0.0};
        auto o2 = quaternion<double> {1.0, 0.0, 0.0, 0.0};
        particle part1(1., 1., p1, o1);
        particle part2(1., 1., p2, o2);
        std::vector<particle> particlelist {part1, part2};

        /* Number of columns in each row of data to be written into HD5 file. It can be 4 or 8.
         * A custom value can also be used, but it has to be explicitly written in the template
         * function write2H5file<rowdim> before compiling. */
        size_t numcols = 8;

        // Choose correct child class of trajectory given the current type of particles
        if (integ.getParticlesBodyType() == "point"){
            traj = std::make_unique<trajectoryPosition>(particlelist.size(), bufferSize);
            //numcols = 4;
        } else {
            traj = std::make_unique<trajectoryPositionOrientation>(particlelist.size(), bufferSize);
            //numcols = 8;
        }

        // Main simulation loop (integration and writing to file)
        for (int tstep=0; tstep < Nsteps; tstep++) {
            integ.integrate(particlelist);
            if (tstep % stride == 0) {
                bufferCounter++;
                traj->sample(tstep, particlelist);
                if (bufferCounter == bufferSize) {
                    bufferCounter = 0;
                    write2H5file(numcols, filename, true);
                    traj->emptyBuffer();
                }
            }
        }

        // Empty remaining data in buffer into H5 file
        if (bufferCounter > 0) {
            write2H5file(numcols, filename, true);
            traj->emptyBuffer();
        }

        // Write to HD5 file
        //if (rowdim == 4) { traj->write2H5file<4>(filename, traj->getData()); }
        //else if (rowdim == 8) { traj->write2H5file<8>(filename, traj->getData()); }
        traj->write2file(filename, traj->getData());
    }


    // Wrapper for traj->write2H5file and traj->writeChunk2H5file
    void simulation::write2H5file(int numcols, std::string filename, bool chunked ) {
        if (chunked) {
            if (numcols == 4) {
                traj->writeChunk2H5file<4>(filename, traj->getData());
            } else if (numcols == 8) {
                traj->writeChunk2H5file<8>(filename, traj->getData());
            } else {
                throw std::range_error("Numcols needs to be 4 or 8. Custom number of columns per row "
                                       "can be used but needs to be explicitly modified in "
                                       "simulation.cpp, write2H5file<numcols> ");
            }
        } else {
            if (numcols == 4) {
                traj->write2H5file<4>(filename, traj->getData());
            } else if (numcols == 8) {
                traj->write2H5file<8>(filename, traj->getData());
            } else {
                throw std::range_error("Numcols needs to be 4 or 8. Custom number of columns per row "
                                       "can be used but needs to be explicitly modified in "
                                       "simulation.cpp, write2H5file<numcols> ");
            }
        }
    }
}