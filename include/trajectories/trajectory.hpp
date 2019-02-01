//
// Created by dibakma on 18.09.18.
//

#pragma once
#include <array>
#include <functional>
#include <vector>
#include <fstream>
#include <iterator>
#include <H5f90i.h>
#include "H5Cpp.h"
#include "particle.hpp"


using namespace std::placeholders;

// Needed to write to HDF5 files.
using namespace H5;

namespace msmrd {
    /**
     * Abstract base class to store full trajectories
     */
    class trajectory {
    protected:
        const std::size_t kB = 1024;
        const std::size_t MB = 1024 * kB;
        bool firstrun = true;
        std::vector<std::vector<double>> trajectoryData;
        std::vector<std::vector<int>> discreteTrajectoryData;
    public:
        unsigned long Nparticles;
        int bufferSize;
        /**
         * @param kB/MB, constant buffer sizes in bytes for writing data
         * @param chunksWritten keeps track of the number of chunks written to file
         * @param Nparticles is the number of particles in the trajectory. In the case of relative
         * sampling of cooridnates, it should correspond to the number of all possible pairs of
         * particles.
         * @bufferSize buffer size for data storage. Exact value if data being dumped into file;
         * otherwise an approximated value is enough.
         * @param trajectoryData buffer to store trajectory data (time, position, and/or other variables
         * like orientation). Defined as a vector of vectors instead of vector of arrays for binding
         * to work.
         * @param discreteTrajectoryData buffer to store the discretized trajectory data (usually in the form
         * of states given by integers). Therefore, defined as a vector of integers.
         */

        trajectory(unsigned long Nparticles, int bufferSize);

        void emptyBuffer();

        // Virtual functions to sample from list of particles, store in trajectoryData, and empty data buffer
        virtual void sample(double time, std::vector<particle> &particleList) = 0;

        virtual void sampleRelative(double time, std::vector<particle> &particleList) = 0;

        virtual void sampleDiscreteTrajectory(double time, std::vector<particle> &particleList) = 0;


        // Functions used by child classes
        std::vector<std::vector<double>> getTrajectoryData() { return trajectoryData; }

        std::vector<std::vector<int>> getDiscreteTrajectoryData() { return discreteTrajectoryData; }

        void write2file(std::string filename, std::vector<std::vector<double>> localdata);

        // Templated functions for writing to H5 file (template useful since datasize needs to be known at runtime)
        template< typename scalar, size_t NUMCOL>
        void write2H5file(std::string filename, std::vector<std::vector<scalar>> localdata);

        template< typename scalar, size_t NUMCOL>
        void writeChunk2H5file(std::string filename, std::vector<std::vector<scalar>> localdata);

    };

    /**
     * Class to store position only trajectories
     */
    class trajectoryPosition: public trajectory {
    public:

        trajectoryPosition(unsigned long Nparticles, int bufferSize);

        void sample(double time, std::vector<particle> &particleList) override;

        void sampleRelative(double time, std::vector<particle> &particleList) override;

        // Empty function to be overwritten by child classes if neccesary
        void sampleDiscreteTrajectory(double time, std::vector<particle> &particleList) override {};

    };


    /**
     * Class to store trajectories with position and orientation (given by a quaternion)
     */
    class trajectoryPositionOrientation : public trajectory {
    public:

        trajectoryPositionOrientation(unsigned long Nparticles, int bufferSize);

        void sample(double time, std::vector<particle> &particleList) override;

        void sampleRelative(double time, std::vector<particle> &particleList) override;

        // Empty function to be overwritten by child classes if neccesary
        void sampleDiscreteTrajectory(double time, std::vector<particle> &particleList) override {};

        void printTime();

    };


    /**
     * Templated trajectory functions for H5 file writing (implementations need to be in header)
     * @param NUMCOL gives the length of elements in each row of data to be written (number of columns)
     */

    // Writes data into HDF5 binary file
    template< typename scalar, size_t NUMCOL>
    void trajectory::write2H5file(std::string filename, std::vector<std::vector<scalar>> localdata) {
        const H5std_string FILE_NAME = filename + ".h5";
        const H5std_string	DATASET_NAME = "msmrd_data";
        int datasize = static_cast<int>(localdata.size());


        // Copies data into fixed size array , datafixed
        double datafixed[datasize][NUMCOL];
        for (int i = 0; i < datasize; i++) {
            for (int j = 0; j < NUMCOL; j++) {
                datafixed[i][j] = 1.0*localdata[i][j];
            }
        }

        // Creates H5 file (overwrites previous existing one)
        H5File file(FILE_NAME, H5F_ACC_TRUNC);

        // Sets shape of data into dataspace
        hsize_t dims[2];               // dataset dimensions
        dims[0] = localdata.size();
        dims[1] = NUMCOL;
        DataSpace dataspace(2, dims);

        // Creates dataset and write data into it
        DataSet dataset = file.createDataSet(DATASET_NAME, PredType::NATIVE_DOUBLE, dataspace);
        dataset.write(datafixed, H5::PredType::NATIVE_DOUBLE);

    };

    // Writes data into HDF5 binary file in chunks of size bufferSize/bufferSize*Nparticles
    template< typename scalar, size_t NUMCOL >
    void trajectory::writeChunk2H5file(std::string filename, std::vector<std::vector<scalar>> localdata) {
        const H5std_string FILE_NAME( filename + ".h5");
        const H5std_string DATASET_NAME( "msmrd_data" );
        hsize_t chunckSize = localdata.size();
        const int RANK = 2;

        H5File file;
        DataSet dataset;
        DataSpace dataspace;
        hsize_t dimsFile[RANK] = {0 ,0};

        // Copies data into fixed size array , datafixed
        double datafixed[chunckSize][NUMCOL];
        for (int i = 0; i < chunckSize; i++) {
            for (int j = 0; j < NUMCOL; j++) {
                datafixed[i][j] = 1.0*localdata[i][j];
            }
        }

        if (firstrun) {
            // Create dataspace with unlimited dimensions
            hsize_t dims[2]  = {chunckSize, NUMCOL};  // dataset dimensions at creation
            hsize_t maxdims[2] = {H5S_UNLIMITED, H5S_UNLIMITED};
            dataspace = DataSpace(RANK , dims, maxdims);

            // Create H5 file. If file exists, it will be overwritten
            file = H5File(FILE_NAME, H5F_ACC_TRUNC);

            // Modify dataset creation properties, i.e. enable chunking.
            DSetCreatPropList cparms;
            hsize_t chunk_dims[2] ={chunckSize, NUMCOL};
            cparms.setChunk( RANK, chunk_dims );

            // Set fill value for the dataset
            double fill_val = 0;
            cparms.setFillValue( PredType::NATIVE_DOUBLE, &fill_val);

            /* Create a new dataset within the file using cparms
            * creation properties.  */
            dataset = file.createDataSet( DATASET_NAME, PredType::NATIVE_DOUBLE, dataspace, cparms);

            firstrun = false;
        }
        else {
            // Open existing dataset
            file = H5File(FILE_NAME, H5F_ACC_RDWR);
            dataset = file.openDataSet(DATASET_NAME);
            dataspace = dataset.getSpace();

            // Get dimensions of current dataset in file
            const int ndims = dataspace.getSimpleExtentDims(dimsFile, NULL);
        }

        // Extend the dataset by a chunk (chunkSize, NUMCOL)
        hsize_t size[2];
        size[0] = dimsFile[0] + chunckSize;
        size[1] = NUMCOL;
        dataset.extend( size );

       // Select a hyperslab.
        DataSpace fspaceChunck = dataset.getSpace ();
        hsize_t offset[2];
        offset[0] = dimsFile[0];
        offset[1] = 0;
        hsize_t dimsChunk[2] = { chunckSize, NUMCOL};            /* data1 dimensions */
        fspaceChunck.selectHyperslab( H5S_SELECT_SET, dimsChunk, offset );

        //Define memory space
        DataSpace mspaceChunk( RANK, dimsChunk );

        // Write the data to the hyperslab.
        dataset.write( datafixed, PredType::NATIVE_DOUBLE, mspaceChunk, fspaceChunck );

    }


}