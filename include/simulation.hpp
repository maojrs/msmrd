//
// Created by dibakma on 22.06.18.
//
#pragma once

#include <array>
#include "particle.hpp"
#include "integrators/integrator.hpp"
#include "trajectories/trajectory.hpp"

namespace msmrd {
    class simulation {
    public:
        std::vector<particle> particleList;
        integrator &integ;
        std::unique_ptr<trajectory> traj;
        int numcols = 4;
        /**
         * @param particleList List of particles to be integrated. In principle it can also take lists of
         * cutom particle types as long as they are childs of the original particle class. However, watch out
         * in the bindings, since custom types might not work.
         * @param integ Integrator to be used for simulation, works for any integrator since they are all
         * childs from abstract class.
         * @param traj smart pointer to trajectory class. The class will be initializaed into one of the
         * child classes of trajectory class.
         * @param numcols Number of columns in each row of data to be written into HD5 file. It can be 4 or 8.
         * A custom value can also be used, but it has to be explicitly written in the template
         * function write2H5file<rowdim> before compiling.
         */


        simulation(integrator &integ);

        void run(int Nsteps, int stride, int bufferSize, const std::string &filename,
                 bool outputTxt, bool outputH5, bool outputChunked);

    private:

        void runNoutputChunks(int Nsteps, int stride, int bufferSize, const std::string &filename, bool chunked);

        void runNoutput(int Nsteps, int stride, int bufferSize, const std::string &filename,
                        bool outputTxt, bool H5output, bool chunked);

        void write2H5file(int numcols, std::string filename, bool chunked); // Wrapper for traj.write2H5file

    };

}