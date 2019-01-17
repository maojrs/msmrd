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
    private:
        const size_t rowdim = 0;
    protected:
        const std::size_t kB = 1024;
        const std::size_t MB = 1024 * kB;
        bool firstrun = true;
        std::vector<std::vector<double>> data;
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
         * @param data buffer to store trajectory data (time, position, and/or other variables
         * like orientation). It will be redefined in child classes to avoid template, so this
         * value remains hidden in child classes)
         */

        trajectory(unsigned long Nparticles, int bufferSize);

        // Virtual functions to sample from list of particles and store in data
        virtual void sample(double time, std::vector<particle> &particleList) = 0;

        virtual void sampleRelative(double time, std::vector<particle> &particleList) = 0;

        virtual void emptyBuffer() = 0;

        std::vector<std::vector<double>> getData() { return data; }

        const size_t getRowDim() { return rowdim; }

        void write2file(std::string filename, std::vector<std::vector<double>> localdata);

        template< size_t ROWDIM>
        void write2H5file(std::string filename, std::vector<std::vector<double>> localdata);

        template< size_t ROWDIM>
        void writeChunk2H5file(std::string filename, std::vector<std::vector<double>> localdata);


        void createExtendibleH5File(std::string filename, hsize_t rowdim);


//        void append(sample_type sample) {
//            data.push_back(sample);
//        };
//
//        sample_type operator[](std::size_t i) const {
//            return data.at(i);
//        };
//
//        sample_type &operator[](std::size_t i) {
//            return data.at(i);
//        };
    };

    /**
     * Class to store position only trajectories
     */
    class trajectoryPosition: public trajectory {
    private:
        const size_t rowdim = 4;
    public:

        trajectoryPosition(unsigned long Nparticles, int bufferSize);

        void sample(double time, std::vector<particle> &particleList) override;

        void sampleRelative(double time, std::vector<particle> &particleList) override;

        void emptyBuffer() override { data.clear(); }

        const size_t getRowDim() { return rowdim; }

        //std::vector<std::array<double, 4>> getData() const { return data; };

    };


    /**
     * Class to store trajectories with position and orientation (given by a quaternion)
     */
    class trajectoryPositionOrientation : public trajectory {
    private:
        const size_t rowdim = 8;
    public:

        trajectoryPositionOrientation(unsigned long Nparticles, int bufferSize);

        void sample(double time, std::vector<particle> &particleList) override;

        void sampleRelative(double time, std::vector<particle> &particleList) override;

        void emptyBuffer() override { data.clear(); }

        const size_t getRowDim() { return rowdim; }

        //std::vector<std::array<double, 8>> getData() const { return data; };

        void printTime();

    };


    /**
     * Templated trajectory functions for H5 file writing (implementations need to be in header)
     * @param ROWDIM gives the length of elements in each row of data to be written
     */

    // Writes data into HDF5 binary file
    template< size_t ROWDIM >
    void trajectory::write2H5file(std::string filename, std::vector<std::vector<double>> localdata) {
        const H5std_string FILE_NAME = filename + ".h5";
        const H5std_string	DATASET_NAME = "msmrd_data";
        int datasize = static_cast<int>(localdata.size());


        // Copies data into fixed size array , datafixed
        double datafixed[datasize][ROWDIM];
        for (int i = 0; i < datasize; i++) {
            for (int j = 0; j < ROWDIM; j++) {
                datafixed[i][j] = 1.0*localdata[i][j];
            }
        }

        // Creates H5 file (overwrites previous existing one)
        H5File file(FILE_NAME, H5F_ACC_TRUNC);

        // Sets shape of data into dataspace
        hsize_t dims[2];               // dataset dimensions
        dims[0] = localdata.size();
        dims[1] = ROWDIM;
        DataSpace dataspace(2, dims);

        // Creates dataset and write data into it
        DataSet dataset = file.createDataSet(DATASET_NAME, PredType::NATIVE_DOUBLE, dataspace);
        dataset.write(datafixed, H5::PredType::NATIVE_DOUBLE);

    };

    // NOTE THIS FUNCTION HASN"T BEEN TESTED, MAY BE INCOMPLETE OR NONFUNCTIONAL
    template< size_t ROWDIM >
    void trajectory::writeChunk2H5file(std::string filename, std::vector<std::vector<double>> localdata) {
        const H5std_string FILE_NAME( filename + ".h5");
        const H5std_string DATASET_NAME( "msmrd_data" );
        hsize_t chunckSize = localdata.size();
        const int RANK = 2;

        H5File file;
        DataSet dataset;
        DataSpace dataspace;
        hsize_t dimsFile[RANK] = {0 ,0};

        // Copies data into fixed size array , datafixed
        double datafixed[chunckSize][ROWDIM];
        for (int i = 0; i < chunckSize; i++) {
            for (int j = 0; j < ROWDIM; j++) {
                datafixed[i][j] = 1.0*localdata[i][j];
            }
        }

        if (firstrun) {
            // Create dataspace with unlimited dimensions
            hsize_t dims[2]  = {chunckSize, ROWDIM};  // dataset dimensions at creation
            hsize_t maxdims[2] = {H5S_UNLIMITED, H5S_UNLIMITED};
            dataspace = DataSpace(RANK , dims, maxdims);

            // Create H5 file. If file exists, it will be overwritten
            file = H5File(FILE_NAME, H5F_ACC_TRUNC);

            // Modify dataset creation properties, i.e. enable chunking.
            DSetCreatPropList cparms;
            hsize_t chunk_dims[2] ={chunckSize, ROWDIM};
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

        // Extend the dataset by a chunk (chunkSize, ROWDIM)
        hsize_t size[2];
        size[0] = dimsFile[0] + chunckSize;
        size[1] = ROWDIM;
        dataset.extend( size );

       // Select a hyperslab.
        DataSpace fspaceChunck = dataset.getSpace ();
        hsize_t offset[2];
        offset[0] = dimsFile[0];
        offset[1] = 0;
        hsize_t dimsChunk[2] = { chunckSize, ROWDIM};            /* data1 dimensions */
        fspaceChunck.selectHyperslab( H5S_SELECT_SET, dimsChunk, offset );

        //Define memory space
        DataSpace mspaceChunk( RANK, dimsChunk );

        // Write the data to the hyperslab.
        dataset.write( datafixed, PredType::NATIVE_DOUBLE, mspaceChunk, fspaceChunck );

    }


}