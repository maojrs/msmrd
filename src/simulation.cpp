//
// Created by dibakma on 22.06.18.
//

//#include <random>
#include "simulation.hpp"
#include <memory>
#include "trajectories/trajectory.hpp"
#include "trajectories/patchyDimer.hpp"



namespace msmrd {
    /*
     * Simulation main class declaration.
     * @param &integ reference to integrator to be used for simulation, works for any integrator since they are all
     * childs from abstract class.
     */
    simulation::simulation(integrator &integ) : integ(integ){};

    // Note the folder where the data is saved should already exist
    void simulation::run(std::vector<particle> &particleList, int Nsteps, int stride, int bufferSize,
                         const std::string &filename, bool outputTxt, bool outputH5, bool outputChunked,
                         std::string trajtype) {

        if (outputTxt && outputChunked){
            throw std::range_error("Output in chunks is not available with txt output. It is recommended to "
                                   "change output to H5 in chunks and turn off txt output.");
        }

        if (!outputH5 && outputChunked) {
            throw std::range_error("Output in chunks is only available with H5 output. It is recommended to "
                                   "change to output with H5 in chunks and turn off txt ouput.");
        }

        // Choose correct child class of trajectory given the current type of particles
        if (trajtype == "patchyDimer") {
            int numSections = 6;
            outputDiscreteTraj = true;
            traj = std::make_unique<patchyDimer>(particleList.size(), bufferSize, 6);
            numcols = 8;
        } else if (trajtype == "position"){
            traj = std::make_unique<trajectoryPosition>(particleList.size(), bufferSize);
            numcols = 4;
        } else if (trajtype == "positionOrientation") {
            traj = std::make_unique<trajectoryPositionOrientation>(particleList.size(), bufferSize);
            numcols = 8;
        } else { // Otherwise use trajectoryPositionOrientation as default class
            traj = std::make_unique<trajectoryPositionOrientation>(particleList.size(), bufferSize);
            numcols = 8;
        }

        /* Main simulation loop. Simulation method depends on output method. If H5 and chunked outputs are chosen
         * the data will be dumped into file everytime the buffer is full and erased from memory. If data is not,
         * chunked, it can be written directyl from memory into a H5 file or a text file, the data is not erased
         * from memory. */
        if (outputH5 && outputChunked) {
            runNoutputChunks(particleList, Nsteps, stride, bufferSize, filename, outputChunked);
        } else {
            runNoutput(particleList, Nsteps, stride, bufferSize, filename, outputTxt, outputH5, outputChunked);
        }

    }


    // Runs simulation while outputing chunked data into H5 file and freeing up memory
    void simulation::runNoutputChunks(std::vector<particle> &particleList, int Nsteps, int stride, int bufferSize,
                                      const std::string &filename, bool chunked){
        int bufferCounter = 0;
        // Main simulation loop (integration and writing to file)
        for (int tstep=0; tstep < Nsteps; tstep++) {
            if (tstep % stride == 0) {
                bufferCounter++;
                traj->sample(tstep, particleList);
                if (outputDiscreteTraj) {
                    traj->sampleDiscreteTrajectory(tstep, particleList);
                }
                if (bufferCounter == bufferSize) {
                    bufferCounter = 0;
                    write2H5file(numcols, filename, true);
                    traj->emptyBuffer();
                }
            }
            integ.integrate(particleList);
        }

        // Empty remaining data in buffer into H5 file
        if (bufferCounter > 0) {
            write2H5file(numcols, filename, true);
            traj->emptyBuffer();
        }
    }

    // Runs simulation, when done outputs data into H5 file , text file or both. Memory is not freed up.
    void simulation::runNoutput(std::vector<particle> &particleList, int Nsteps, int stride, int bufferSize,
                                const std::string &filename, bool outputTxt, bool outputH5, bool chunked){
        // Main simulation loop (integration and writing to file)
        for (int tstep=0; tstep < Nsteps; tstep++) {
            if (tstep % stride == 0) {
                traj->sample(tstep, particleList);
            }
            integ.integrate(particleList);
        }
        // Writes into H5 file
        if (outputH5){
            write2H5file(numcols, filename, chunked);
        }
        // writes into normal textfile
        if (outputTxt) {
            traj->write2file(filename, traj->getTrajectoryData());
        }
    }


    // Wrapper for traj->write2H5file and traj->writeChunk2H5file
    void simulation::write2H5file(int numcols, std::string filename, bool chunked ) {
        if (chunked) {
            // Write discrete trajectory chunked
//            if (outputDiscreteTraj) {
//                traj->writeChunk2H5file<int, 1>(filename + "_discrete", traj->getDiscreteTrajectoryData());
//            }
            // Write the continuous trajectory chunked
            if (numcols == 4) {
                traj->writeChunk2H5file<double, 4>(filename, traj->getTrajectoryData());
            } else if (numcols == 8) {
                traj->writeChunk2H5file<double, 8>(filename, traj->getTrajectoryData());
            } else {
                throw std::range_error("Numcols needs to be 4 or 8. Custom number of columns per row "
                                       "can be used but needs to be explicitly modified in "
                                       "simulation.cpp, write2H5file<numcols> ");
            }
        } else {
            // Write discrete trajectory
//            if (outputDiscreteTraj) {
//                traj->write2H5file<int, 1>(filename + "_discrete", traj->getDiscreteTrajectoryData());
//            }
            // Write the continuous trajectory
            if (numcols == 4) {
                traj->write2H5file<double, 4>(filename, traj->getTrajectoryData());
            } else if (numcols == 8) {
                traj->write2H5file<double, 8>(filename, traj->getTrajectoryData());
            } else {
                throw std::range_error("Numcols needs to be 4 or 8. Custom number of columns per row "
                                       "can be used but needs to be explicitly modified in "
                                       "simulation.cpp, write2H5file<numcols> ");
            }
        }
    }
}