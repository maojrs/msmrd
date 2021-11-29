//
// Created by dibakma on 22.06.18.
//
#pragma once

#include <array>
#include <memory>
#include "particle.hpp"
#include "integrators/integrator.hpp"
#include "trajectories/trajectory.hpp"
#include "trajectories/trajectoryEnergyTemperature.hpp"
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
        std::unique_ptr<trajectory> trajEnergyTemp;
        std::vector<int> distinguishedTypes{1};
        int numcols = 4;
        bool outputDiscreteTraj = false;
        bool outputEnergyTemperature = false;
        int equilibrationSteps = 0;
        /**
         * @param integ Integrator to be used for simulation, works for any integrator since they are all
         * childs from abstract class.
         * @param traj smart pointer to trajectory class. The class will be initializaed into one of the
         * child classes of trajectory class.
         * @param trajEnergyTemp smart pointer to trjaectory class for energy and temperature. Only available for
         * Langevin simulations and only used if outputEnergyTemperature = true.
         * @param distinguishedTypes Set distinguished types for distinguished trajectories (1 by default).
         * This means distinguished trajectories will only sample particles of type 1
         * @param numcols Number of columns in each row of data to be written into HD5 file. It can be 4 or 8.
         * A custom value can also be used, but it has to be explicitly written in the template
         * function write2H5file<rowdim> before compiling.
         * @param outputDiscreteTraj if true, outputs discrete trajectory. Only available for certain
         * trajectory classes.
         * @param outputEnergyTemperature only valid for Langevin trajectories, if true outputs energy and
         * instant temperature of system at every time step in a different file.
         */


        simulation(integrator &integ);

        void run(std::vector<particle> &particleList, int Nsteps, int stride, int bufferSize,
                 const std::string &filename, bool outputTxt, bool outputH5, bool outputChunked,
                 std::string trajtype);

        void setEquilibrationSteps(int eqSteps) {
            equilibrationSteps = eqSteps;
        }

        void setOutputEnergyTemperature(bool value) {
            outputEnergyTemperature = value;
        }

        void setDistinguishedTypes(std::vector<int> distTypes) {
            distinguishedTypes.resize(distTypes.size());
            distinguishedTypes = distTypes;
        }

    private:

        void runNoutputChunks(std::vector<particle> &particleList, int Nsteps, int stride, int bufferSize,
                              const std::string &filename);

        void runNoutput(std::vector<particle> &particleList, int Nsteps, int stride, int bufferSize,
                        const std::string &filename, bool outputTxt, bool H5output);

        void runEquilibration(std::vector<particle> &particleList);

        void createChunkedH5files(std::string filename);

        void write2H5file(std::string filename, bool chunked); // Wrapper for traj.write2H5file

    };

}