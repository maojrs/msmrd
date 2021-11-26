//
// Created by dibakma on 22.06.18.
//

//#include <random>
#include "simulation.hpp"



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
            throw std::invalid_argument("Output in chunks is not available with txt output. It is recommended to "
                                        "change output to H5 in chunks and turn off txt output.");
        }

        if (!outputH5 && outputChunked) {
            throw std::invalid_argument("Output in chunks is only available with H5 output. It is recommended to "
                                        "change to output with H5 in chunks and turn off txt ouput.");
        }

        /* Set distinguished typer for distinguished trajectories (1 by default). This means distinguished trajectories
         * will only sample particles of type 1 */
        std::vector<int> distinguishedTypes{1};

        // Choose correct child class of trajectory given the current type of particles
        if (trajtype == "patchyDimer") {
            outputDiscreteTraj = false;
            traj = std::make_unique<patchyDimerTrajectory>(particleList.size(), bufferSize);
            numcols = 9; //(time, positionx3, orientationx4, state)
        } else if (trajtype == "patchyDimer2"){
            outputDiscreteTraj = false;
            traj = std::make_unique<patchyDimerTrajectory2>(particleList.size(), bufferSize);
            numcols = 9; //(time, positionx3, orientationx4, state)
        } else if (trajtype == "patchyProtein") {
            outputDiscreteTraj = false;
            traj = std::make_unique<patchyProteinTrajectory>(particleList.size(), bufferSize);
            numcols = 9; //(time, positionx3, orientationx4, state)
        } else if (trajtype == "patchyProtein2") {
            outputDiscreteTraj = false;
            traj = std::make_unique<patchyProteinTrajectory2>(particleList.size(), bufferSize);
            numcols = 9; //(time, positionx3, orientationx4, state)
        } else if (trajtype == "MAPK"){
            outputDiscreteTraj = false;
            traj = std::make_unique<MAPKtrajectory>(particleList.size(), bufferSize);
            numcols = 10; //(time, positionx3, orientationx4, state, type)
        } else if (trajtype == "moriZwanzig") {
            outputDiscreteTraj = false;
            traj = std::make_unique<trajectoryMoriZwanzig>(particleList.size(), bufferSize,
                                                           distinguishedTypes);
            numcols = 8; //(time, positionx3, type, raux(x3))
        } else if (trajtype == "moriZwanzigVelocity"){
                outputDiscreteTraj = false;
                traj = std::make_unique<trajectoryMoriZwanzigVelocity>(particleList.size(), bufferSize,
                                                               distinguishedTypes);
                numcols = 11; //(time, positionx3, velocityx3, type, raux(x3))
        } else if (trajtype == "position"){
            traj = std::make_unique<trajectoryPosition>(particleList.size(), bufferSize);
            numcols = 4; //(time, positionx3)
        }  else if (trajtype == "positionType"){
            traj = std::make_unique<trajectoryPositionType>(particleList.size(), bufferSize);
            numcols = 5; //(time, positionx3, type)
        } else if (trajtype == "positionDistinguished"){
            traj = std::make_unique<trajectoryPositionDistinguished>(particleList.size(), bufferSize,
                                                                     distinguishedTypes);
            numcols = 5; //(time, positionx3, type)
        } else if (trajtype == "positionVelocity"){
            traj = std::make_unique<trajectoryPositionVelocity>(particleList.size(), bufferSize);
            numcols = 7; //(time, positionx3, velocityx3)
        } else if (trajtype == "positionVelocityType"){
            traj = std::make_unique<trajectoryPositionVelocityType>(particleList.size(), bufferSize);
            numcols = 8; //(time, positionx3, velocityx3, type)
        } else if (trajtype == "positionVelocityDistinguished"){
            traj = std::make_unique<trajectoryPositionVelocityDistinguished>(particleList.size(), bufferSize,
                                                                    distinguishedTypes);
            numcols = 8; //(time, positionx3, velocityx3, type)
        } else if (trajtype == "positionOrientation") {
            traj = std::make_unique<trajectoryPositionOrientation>(particleList.size(), bufferSize);
            numcols = 8; //(time, positionx3, orientationx4)
        } else if (trajtype == "positionOrientationType") {
            traj = std::make_unique<trajectoryPositionOrientationType>(particleList.size(), bufferSize);
            numcols = 9; //(time, positionx3, orientationx4)
        } else if (trajtype == "positionOrientationState") {
            traj = std::make_unique<trajectoryPositionOrientationState>(particleList.size(), bufferSize);
            numcols = 9; //(time, positionx3, orientationx4, state)
        } else if (trajtype == "positionOrientationStateType") {
            traj = std::make_unique<trajectoryPositionOrientationStateType>(particleList.size(), bufferSize);
            numcols = 10; //(time, positionx3, orientationx4, state, type)
        } else { // Otherwise use trajectoryPositionOrientationStateType as default class
            traj = std::make_unique<trajectoryPositionOrientationStateType>(particleList.size(), bufferSize);
            numcols = 10; ////(time, positionx3, orientationx4, state, Type)
        }

        // Set boundary in trajectory class
        if (integ.isBoundaryActive()) {
            traj->setBoundary(integ.getBoundary());
        }

        // Set energy temperature trajectory
        if (outputEnergyTemperature) {
            trajEnergyTemp = std::make_unique<trajectoryEnergyTemperature>(bufferSize);
            if (integ.isExternalPotentialActive()) {
                trajEnergyTemp->setExternalPotential(integ.getExternalPotential());
            }
            if (integ.isPairPotentialActive()) {
                trajEnergyTemp->setPairPotential(integ.getPairPotential());
            }
        }

        /* Equilibration step (run system for equilibrationStates timsteps) in case it is needed
         * to equilibrate the sytsem before the main simulation. The default value is zero
         * equilibration steps, so one needs to set the equilibrationSteps first, see setEquilibrationSteps. */
        if (equilibrationSteps != 0) {
            runEquilibration(particleList);
        }

        /* Main simulation loop. Simulation method depends on output method. If H5 and chunked outputs are chosen
         * the data will be dumped into file everytime the buffer is full and erased from memory. If data is not,
         * chunked, it can be written directyl from memory into a H5 file or a text file, the data is not erased
         * from memory. */
        if (outputH5 && outputChunked) {
            runNoutputChunks(particleList, Nsteps, stride, bufferSize, filename);
        } else {
            runNoutput(particleList, Nsteps, stride, bufferSize, filename, outputTxt, outputH5);
        }

    }


    // Runs simulation while outputing chunked data into H5 file and freeing up memory
    void simulation::runNoutputChunks(std::vector<particle> &particleList, int Nsteps, int stride, int bufferSize,
                                      const std::string &filename){
        int bufferCounter = 0;
        bool chunked = true;
        // Main simulation loop (integration and writing to file)
        for (int tstep=0; tstep < Nsteps; tstep++) {
            if (tstep % stride == 0) {
                bufferCounter++;
                traj->sample(integ.clock, particleList);
                if (outputDiscreteTraj) {
                    traj->sampleDiscreteTrajectory(integ.clock, particleList);
                }
                if (outputEnergyTemperature) {
                    trajEnergyTemp->sample(integ.clock, particleList);
                }
                if (tstep == 0) {
                    createChunkedH5files(filename);
                }

                //Write to file once buffer is full
                if (bufferCounter * particleList.size() >= bufferSize) {
                    bufferCounter = 0;
                    write2H5file(filename, chunked);
                    traj->emptyBuffer();
                    if (outputEnergyTemperature) {
                        trajEnergyTemp->emptyBuffer();
                    }
                }
            }
            integ.integrate(particleList);
        }

        // Empty remaining data in buffer into H5 file
        if (bufferCounter > 0) {
            write2H5file(filename, chunked);
            traj->emptyBuffer();
            if (outputEnergyTemperature) {
                trajEnergyTemp->emptyBuffer();
            }
        }
    }

    // Runs simulation, when done outputs data into H5 file , text file or both. Memory is not freed up.
    void simulation::runNoutput(std::vector<particle> &particleList, int Nsteps, int stride, int bufferSize,
                                const std::string &filename, bool outputTxt, bool outputH5){
        bool chunked = false;
        // Main simulation loop (integration and writing to file)
        for (int tstep=0; tstep < Nsteps; tstep++) {
            if (tstep % stride == 0) {
                traj->sample(integ.clock, particleList);
            }
            integ.integrate(particleList);
        }
        // Writes into H5 file
        if (outputH5){
            write2H5file(filename, chunked);
        }
        // writes into normal textfile
        if (outputTxt) {
            traj->write2file<double>(filename, traj->getTrajectoryData());
            if (outputDiscreteTraj) {
                traj->write2file<int>(filename + "_discrete", traj->getDiscreteTrajectoryData());
            }
        }
    }

    // Runs integrator for equilibrationSteps, timesteps to equilibarte the system
    void simulation::runEquilibration(std::vector<particle> &particleList) {
        // Main simulation loop (integration and writing to file)
        for (int tstep=0; tstep < equilibrationSteps; tstep++) {
            integ.integrate(particleList);
        }
    }

    // Wrapper for creating chunkedH5files
    void simulation::createChunkedH5files(std::string filename) {
        // Create H5 files to refill discrete trajectory by chunks
        if (outputDiscreteTraj) {
            traj->createChunkedH5file<int, 1>(filename + "_discrete", "msmrd_discrete_data",
                                              traj->getDiscreteTrajectoryData());
        }
        // Create H5 files to refill energy temperature trajectory by chunks
        if (outputEnergyTemperature) {
            trajEnergyTemp->createChunkedH5file<double, 3>(filename + "_energytemp", "msmrd_energytemp_data",
                                              trajEnergyTemp->getTrajectoryData());
        }
        // Create H5 files to refill trajectory by chunks
        if (numcols == 4) {
            traj->createChunkedH5file<double, 4>(filename, "msmrd_data", traj->getTrajectoryData());
        } else if (numcols == 5) {
            traj->createChunkedH5file<double, 5>(filename, "msmrd_data", traj->getTrajectoryData());
        } else if (numcols == 6) {
            traj->createChunkedH5file<double, 6>(filename, "msmrd_data", traj->getTrajectoryData());
        } else if (numcols == 7) {
            traj->createChunkedH5file<double, 7>(filename, "msmrd_data", traj->getTrajectoryData());
        } else if (numcols == 8) {
            traj->createChunkedH5file<double, 8>(filename, "msmrd_data", traj->getTrajectoryData());
        } else if (numcols == 9) {
            traj->createChunkedH5file<double, 9>(filename, "msmrd_data", traj->getTrajectoryData());
        } else if (numcols == 10) {
            traj->createChunkedH5file<double, 10>(filename, "msmrd_data", traj->getTrajectoryData());
        } else if (numcols == 11) {
            traj->createChunkedH5file<double, 11>(filename, "msmrd_data", traj->getTrajectoryData());
        } else {
            throw std::invalid_argument("Numcols needs to be an interger in the interval [4,11]. Custom number of columns per row "
                                        "can be used but needs to be explicitly modified in "
                                        "simulation.cpp, write2H5file<numcols> ");
        }
    }

    // Wrapper for traj->write2H5file and traj->writeChunk2H5file
    void simulation::write2H5file(std::string filename, bool chunked ) {
        if (chunked) {
            // Write discrete trajectory chunked
            if (outputDiscreteTraj) {
                traj->writeChunk2H5file<int, 1>(filename + "_discrete", "msmrd_discrete_data", traj->getDiscreteTrajectoryData());
            }
            // Write energy temperature chunked
            if (outputEnergyTemperature) {
                trajEnergyTemp->writeChunk2H5file<double, 3>(filename + "_energytemp", "msmrd_energytemp_data", trajEnergyTemp->getTrajectoryData());
            }
            // Write the continuous trajectory chunked
            if (numcols == 4) {
                traj->writeChunk2H5file<double, 4>(filename, "msmrd_data", traj->getTrajectoryData());
            } else if (numcols == 5) {
                traj->writeChunk2H5file<double, 5>(filename, "msmrd_data", traj->getTrajectoryData());
            } else if (numcols == 6) {
                traj->writeChunk2H5file<double, 6>(filename, "msmrd_data", traj->getTrajectoryData());
            } else if (numcols == 7) {
                traj->writeChunk2H5file<double, 7>(filename, "msmrd_data", traj->getTrajectoryData());
            } else if (numcols == 8) {
                traj->writeChunk2H5file<double, 8>(filename, "msmrd_data", traj->getTrajectoryData());
            } else if (numcols == 9) {
                traj->writeChunk2H5file<double, 9>(filename, "msmrd_data", traj->getTrajectoryData());
            } else if (numcols == 10) {
                traj->writeChunk2H5file<double, 10>(filename, "msmrd_data", traj->getTrajectoryData());
            } else if (numcols == 11) {
                traj->writeChunk2H5file<double, 11>(filename, "msmrd_data", traj->getTrajectoryData());
            } else {
                throw std::invalid_argument("Numcols needs to be an integer in the interval [4,11]. Custom number of columns per row "
                                            "can be used but needs to be explicitly modified in "
                                            "simulation.cpp, write2H5file<numcols> ");
            }
        } else {
            // Write discrete trajectory
            if (outputDiscreteTraj) {
                traj->write2H5file<int, 1>(filename + "_discrete", "msmrd_discrete_data", traj->getDiscreteTrajectoryData());
            }
            // Write energy temperature
            if (outputEnergyTemperature) {
                trajEnergyTemp->write2H5file<double, 3>(filename + "_energytemp", "msmrd_energytemp_data", trajEnergyTemp->getTrajectoryData());
            }
            // Write the continuous trajectory
            if (numcols == 4) {
                traj->write2H5file<double, 4>(filename, "msmrd_data", traj->getTrajectoryData());
            } else if (numcols == 5) {
                traj->write2H5file<double, 5>(filename, "msmrd_data", traj->getTrajectoryData());
            } else if (numcols == 6) {
                traj->write2H5file<double, 6>(filename, "msmrd_data", traj->getTrajectoryData());
            } else if (numcols == 7) {
                traj->write2H5file<double, 7>(filename, "msmrd_data", traj->getTrajectoryData());
            } else if (numcols == 8) {
                traj->write2H5file<double, 8>(filename, "msmrd_data", traj->getTrajectoryData());
            } else if (numcols == 9) {
                traj->write2H5file<double, 9>(filename, "msmrd_data", traj->getTrajectoryData());
            } else if (numcols == 10) {
                traj->write2H5file<double, 10>(filename, "msmrd_data", traj->getTrajectoryData());
            } else if (numcols == 11) {
                traj->write2H5file<double, 11>(filename, "msmrd_data", traj->getTrajectoryData());
            } else {
                throw std::invalid_argument("Numcols needs to be an integer in the interval [4,11]. Custom number of columns per row "
                                            "can be used but needs to be explicitly modified in "
                                            "simulation.cpp, write2H5file<numcols> ");
            }
        }
    }
}
