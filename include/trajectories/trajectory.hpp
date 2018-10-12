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
    public:
        int Nparticles;
        int bufferSize;
        std::vector<std::array<double, 0>> data;
        /**
         * @param kB/MB, constant buffer sizes in bytes for writing data
         * @param Nparticles is the number of particles in the trajectory. In the case of relative
         * sampling of cooridnates, it should correspond to the number of all possible pairs of
         * particles.
         * @bufferSize buffer size for data storage. Exact value if data being dumped into file;
         * otherwise an approximated value is enough.
         * @param data buffer to store trajectory data (time, position, and/or other variables
         * like orientation). It will be redefined in child classes to avoid template, so this
         * value remains hidden in child classes)
         */

        trajectory(int Nparticles, int bufferSize);

        // Virtual functions to sample from list of particles and store in data
        virtual void sample(double time, std::vector<particle> &particleList) = 0;

        virtual void sampleRelative(double time, std::vector<particle> &particleList) = 0;

        virtual void emptyBuffer() = 0;

        template< size_t ROWDIM>
        void write2file(std::string filename, std::vector<std::array<double, ROWDIM>> localdata);

        template< size_t ROWDIM>
        void write2H5file(std::string filename, std::vector<std::array<double, ROWDIM>> localdata);

        template< size_t ROWDIM>
        void write2ExtendibleH5file(std::string filename, std::vector<std::array<double, ROWDIM>> localdata);


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
    class trajectoryPosition : public trajectory {
    private:
        std::vector<std::array<double, 4>> data;
    public:

        trajectoryPosition(int Nparticles, int bufferSize);

        void sample(double time, std::vector<particle> &particleList) override;

        void sampleRelative(double time, std::vector<particle> &particleList) override;

        void emptyBuffer() override { data.clear(); }

        std::vector<std::array<double, 4>> getData() const { return data; };

    };


    /**
     * Class to store trajectories with position and orientation (given by a quaternion)
     */
    class trajectoryPositionOrientation : public trajectory {
    private:
        std::vector<std::array<double, 8>> data;
    public:

        trajectoryPositionOrientation(int Nparticles, int bufferSize);

        void sample(double time, std::vector<particle> &particleList) override;

        void sampleRelative(double time, std::vector<particle> &particleList) override;

        void emptyBuffer() override { data.clear(); }

        std::vector<std::array<double, 8>> getData() const { return data; };

        void printTime();
//        std::function<void(double, std::vector<particle>&)> f = std::bind(&trajectoryPositionOrientation::sample, this, _1, _2);
//        std::function<void(double, std::vector<particle>&)> get_sampler() {
//            return f;
//        }
    };


    /**
     * Templated trajcetory functions for file writing (implementations need to be in header)
     * @param ROWDIM gives the length of elements in each row of data to be written
     */

    // Writes data into normal text file
    template< size_t ROWDIM >
    void trajectory::write2file(std::string filename, std::vector<std::array<double, ROWDIM>> localdata) {
        std::ofstream outputfile(filename + ".txt");
        std::ostream_iterator<double> output_iterator(outputfile, " ");
        for (auto const &value: localdata) {
            std::copy(value.begin(), value.end(), output_iterator);
            outputfile << std::endl;
        }
        outputfile.close();
    };

    // Writes data into HDF5 binary file
    template< size_t ROWDIM >
    void trajectory::write2H5file(std::string filename, std::vector<std::array<double, ROWDIM>> localdata) {
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

        /* Empties data buffer once it has been written to file (note the capacity
         * from the initial data.reserve() remains the same) */
        //data.clear();

    };

    template< size_t ROWDIM >
    void trajectory::write2ExtendibleH5file(std::string filename, std::vector<std::array<double, ROWDIM>> localdata) {
        const H5std_string FILE_NAME( filename + ".h5");
        const H5std_string DATASET_NAME( "msmrd_data" );
        //int datasize = static_cast<int>(localdata.size());
        hsize_t chunckSize = localdata.size();
        const int RANK = 2;

        // Copies data into fixed size array , datafixed
        double datafixed[chunckSize][ROWDIM];
        for (int i = 0; i < chunckSize; i++) {
            for (int j = 0; j < ROWDIM; j++) {
                datafixed[i][j] = 1.0*localdata[i][j];
            }
        }

        // Create dataspace with unlimited dimensions
        hsize_t dims[2]  = {chunckSize, ROWDIM};  // dataset dimensions at creation
        hsize_t maxdims[2] = {H5S_UNLIMITED, H5S_UNLIMITED};
        DataSpace mspace1(RANK , dims, maxdims);

        // Create H5 file. If file exists, it will be overwritten
        H5File file(FILE_NAME, H5F_ACC_TRUNC );

        // Modify dataset creation properties, i.e. enable chunking.
        DSetCreatPropList cparms;
        hsize_t chunk_dims[2] ={chunckSize, ROWDIM};
        cparms.setChunk( RANK, chunk_dims );

        // Set fill value for the dataset
        double fill_val = 0;
        cparms.setFillValue( PredType::NATIVE_DOUBLE, &fill_val);

        /* Create a new dataset within the file using cparms
        * creation properties.  */
        DataSet dataset = file.createDataSet( DATASET_NAME, PredType::NATIVE_DOUBLE, mspace1, cparms);

        // Extend the dataset. This call assures that dataset is at least 3 x 3.
        hsize_t size[2];
        size[0] = chunckSize;
        size[1] = ROWDIM;
        dataset.extend( size );

       // Select a hyperslab.
        DataSpace fspace1 = dataset.getSpace ();
        hsize_t     offset[2];
        offset[0] = 0;
        offset[1] = 0;
        hsize_t dims1[2] = { chunckSize, ROWDIM};            /* data1 dimensions */
        fspace1.selectHyperslab( H5S_SELECT_SET, dims1, offset );

       // Write the data to the hyperslab.
       dataset.write( datafixed, PredType::NATIVE_DOUBLE, mspace1, fspace1 );
    }

} //namespace msmrd