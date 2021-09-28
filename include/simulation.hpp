//
// Created by dibakma on 22.06.18.
//
#pragma once

#include <array>
#include <memory>
#include "particle.hpp"
#include "integrators/integrator.hpp"
#include "trajectories/trajectory.hpp"
#include "trajectories/trajectoryPosition.hpp"
#include "trajectories/trajectoryPositionOrientation.hpp"
#include "trajectories/trajectoryMoriZwanzig.hpp"
#include "trajectories/discrete/patchyDimerTrajectory.hpp"
#include "trajectories/discrete/patchyProteinTrajectory.hpp"
#include "trajectories/discrete/MAPKtrajectory.hpp"


namespace msmrd {
    class simulation {
    public:
        integrator &integ;
        std::unique_ptr<trajectory> traj;
        int numcols = 4;
        bool outputDiscreteTraj = false;
        int equilibrationSteps = 0;
        /**
         * @param integ Integrator to be used for simulation, works for any integrator since they are all
         * childs from abstract class.
         * @param traj smart pointer to trajectory class. The class will be initializaed into one of the
         * child classes of trajectory class.
         * @param numcols Number of columns in each row of data to be written into HD5 file. It can be 4 or 8.
         * A custom value can also be used, but it has to be explicitly written in the template
         * function write2H5file<rowdim> before compiling.
         * @param outputDiscreteTraj if true, outputs discrete trajectory. Only available for certain
         * trajectory classes.
         */


        simulation(integrator &integ);

        void run(std::vector<particle> &particleList, int Nsteps, int stride, int bufferSize,
                 const std::string &filename, bool outputTxt, bool outputH5, bool outputChunked,
                 std::string trajtype);

        void setEquilibrationSteps(int eqSteps) {
            equilibrationSteps = eqSteps;
        }

    private:

        void runNoutputChunks(std::vector<particle> &particleList, int Nsteps, int stride, int bufferSize,
                              const std::string &filename, bool chunked);

        void runNoutput(std::vector<particle> &particleList, int Nsteps, int stride, int bufferSize,
                        const std::string &filename, bool outputTxt, bool H5output, bool chunked);

        void runEquilibration(std::vector<particle> &particleList);

        void createChunkedH5files(std::string filename);

        void write2H5file(std::string filename, bool chunked); // Wrapper for traj.write2H5file

    };

}