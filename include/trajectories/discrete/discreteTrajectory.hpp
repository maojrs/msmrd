//
// Created by maojrs on 9/2/19.
//

#pragma once
#include <memory>
#include <cmath>
#include "trajectories/trajectoryPositionOrientation.hpp"
#include "discretizations/positionOrientationPartition.hpp"
#include "tools.hpp"
#include "H5Cpp.h"

namespace msmrd {
    /**
     * Base class to create discrete trajectories. The main difference with the general class
     * trajectoryPositionOrientation is that this class can sample the discrete trajectory that is
     * specific for the application. In short words, it chooses how to discretize the full
     * trajectory of two particles into a discretized trajectory to be analyzed and extracted into
     * a Markov state model. In general the discretizations will follow the core MSM approach. It
     * also implements functionality to discretize trajectories directly loaded from a python array.
     *
     * @tparam numBoundStates number of bound states in current trajectory discretization. Note
     * the max number of bound states is calculated from this value (8 bound states = 10 max
     * number bound states; 74 states = 100 max number of bound states). This is useful to start
     * the indexing of transition states from a "nice number".
     */

    template<int numBoundStates>
    class discreteTrajectory : public trajectoryPositionOrientation {
    protected:
        std::unique_ptr<positionOrientationPartition> positionOrientationPart;
        std::array< std::tuple<vec3<double>, quaternion<double>>, numBoundStates> boundStates{};
        int maxNumberBoundStates = int (std::pow(10, std::ceil(std::log10(numBoundStates))));
        double rLowerBound = 1.25;
        double rUpperBound = 2.25;
        double tolerancePosition = 0.12;
        double toleranceOrientation = 0.12*2*M_PI;
        int prevsample = 0;
    public:
        /*
         * @positionOrientationPart full six dimensional partition of phase space of relative position and orientation.
         * This is required to sample the discrete trajectory in the transition regions given by this discretization.
         * IMPORTANT: The MSMRD integrator must use the same partition. Also this pointer should be set in constructor
         * of child class.
         * @boundStates vector of tuples. Each tuple contains a vector and a quaternion indicating each one of the
         * bound states. The vector corresponds to the relative position (pos2-pos1) in the frame of reference of
         * particle 1. (smaller index) between particle 1 and 2. The fixed frame of reference also assumes particle 1
         * is in its default initial orientation. The quaternion corresponds to the relative orientation between the
         * two particles, also measured from the fixed frame of reference of particle 1. These bound states are
         * calculated by the setBoundStates function.
         * @param maxNumberBoundStates maximum number of bound states supported. It is used to determine how to
         * count (index) the transition states. The state maxNumberBoundStates + 1 will correspond not to a bound state
         * but to the first transition state. This parameter has to be consistent with the one used
         * by the msmrd integrator and the msmrdMarkovModel.
         * @param rLowerBound any relative distance smaller than this value will no longer assign states using
         * the positionOrientationPartition. Instead in the region bounded by r<rLowerBound, either a bound state
         * is assigned or the Core MSM approach will determine the state opf the discrete trajectory, i.e. the previous
         * state will be sampled in the discrete trajectory until a new bound state is reached or r>=rLowerBound.
         * Note rLowerBound < rUpperBound (positionOrientationPart->relativeDistanceCutOff).
         * @param rUpperBound is the upper relative distance limit for the positionOrientationPartition
         * (positionOrientationPart->relativeDistanceCutOff). Any relative distance above will yield the unbound
         * state 0, instead of a transition state in the region rLowerBound <= r < rUpperBound.
         * IMPORTANT: The MSMRD integrator must have the this same value as the relativeDistanceCutOff.
         * @param tolerancePosition is the maximum acceptable difference between the relative position and the
         * calculated relative position of a metstable region to still be considered part of a bound state.
         * @param toleranceOrientation is the maximum acceptable angle-distance difference between the relative
         * orientation and the relative orientation calculated of a given metastable region to still be considerer
         * part of a bound state.
         * @param prevsample keeps calue of previous sample when sampling discrete trajectory, useful for
         * CoreMSM approach. The CoreMSM approach chooses how to discretize the region r<rLowerBound that
         * is not a bound state. CoreMSM uses the value of the previous known bound or transition state until a new
         * bound or transition state is reached.
         */

        discreteTrajectory(unsigned long Nparticles, int bufferSize);

        discreteTrajectory(unsigned long Nparticles, int bufferSize, double rLowerBound, double rUpperBound);

        void sampleDiscreteTrajectory(double time, std::vector<particle> &particleList) override;

        virtual int sampleDiscreteState(particle part1, particle part2);

        int getBoundState(vec3<double> relativePosition, quaternion<double> relativeOrientation);



        // Virtual function that needs to be overriden by child classes

        virtual void setBoundStates() = 0;


        // Functions are mostly only used when interacting with python or by pybind.

        std::vector<double> discretizeTrajectory(std::vector<std::vector<double>> trajectory);

        // Load H5 directly and discretizes it
        std::vector<double> discretizeTrajectoryH5(std::string filename);

        // Setter functions so child classes can modify default values of parameters

        void setRadialBounds(double rlower, double rupper);

        void setTolerances(double positionTolerance, double orientationTolerance);

    };


    /*
     * Template implementations need to remain in header file, so the whole implementation is below
     */

    /*
     * Chooses how to discretize the full trajectory (trajectoryData) into a discretized trajectory
     * to be analyzed and extracted into a Markov state model. This specific discretization will follow
     * the core MSM approach
     */
    template<int numBoundStates>
    discreteTrajectory<numBoundStates>::discreteTrajectory(unsigned long Nparticles, int bufferSize) :
            trajectoryPositionOrientation(Nparticles, bufferSize) {
        discreteTrajectoryData.reserve(bufferSize);
    };

    template<int numBoundStates>
    discreteTrajectory<numBoundStates>::discreteTrajectory(unsigned long Nparticles, int bufferSize,
                                                           double rLowerBound, double rUpperBound) :
            trajectoryPositionOrientation(Nparticles, bufferSize), rLowerBound(rLowerBound), rUpperBound(rUpperBound) {
        discreteTrajectoryData.reserve(bufferSize);
    };

    /* Samples discrete trajectory directly from the particle list. This is done for every timestep and saved into
     * the class data directly, so this function is better suited to obtain the discrete trajectory at the same
     * time as it is being computed. It is assumed only the first two particles are relevant. */
    template<int numBoundStates>
    void discreteTrajectory<numBoundStates>::sampleDiscreteTrajectory(double time,
                                                                      std::vector<particle> &particleList) {
        // Sample discrete state (use "this->" to make sure it calls the virtual overriden method in child classes).
        int sample = this->sampleDiscreteState(particleList[0], particleList[1]);

        // Apply coreMSM approach
        if (sample == -1) {
            sample = 1*prevsample;
        }

        // Save previous value and push into trajectory
        prevsample = 1*sample;
        discreteTrajectoryData.push_back(std::vector<int>{sample});
    };


    /* Main function to sample the discrete state of two particles. It returns the corresponding
     * bound state, transition state or unbound state (0). In the bound region (r< rLowerBound), it can
     * also return -1 when not in any bound state. In this case, one would normally apply the coreMSM
     * approach and choose the previous value. However, this is done directly on sampleDiscreteTrajectory or
     * in discretizeTrajectoryH5 and discretizeTrajectory if discretizing directly a python array. This function
     * is set a svirtual since it is likely the one that needs to be modified in child classes. */
    template<int numBoundStates>
    int discreteTrajectory<numBoundStates>::sampleDiscreteState(particle part1, particle part2) {
        // Initialize sample with value zero (unbound state)
        int discreteState = 0;

        /* Calculate relative position taking into account periodic boundary measured
         * from i to j (gets you from i to j). */
        vec3<double> relativePosition = calculateRelativePosition(part1.position, part2.position);

        // Rotate relative position to match the reference orientation of particle 1. (VERY IMPORTANT)
        relativePosition = msmrdtools::rotateVec(relativePosition, part1.orientation.conj());
        quaternion<double> quatReference = {1,0,0,0}; // we can then define reference quaternion as identity.

        // Calculate relative orientation (w/respect to particle 1)
        quaternion<double> relativeOrientation;
        //relativeOrientation = part1.orientation.conj() * part2.orientation;
        relativeOrientation =  part2.orientation * part1.orientation.conj();


        // Extract current state, save into sample and return sample
        int secNum;
        if (relativePosition.norm() < rLowerBound) {
            // Returns discrete state or -1 if it is no in any bound state
            discreteState = getBoundState(relativePosition, relativeOrientation);
        }
        // Returns a transitions state if it is in the transition region
        else if (relativePosition.norm() < positionOrientationPart->relativeDistanceCutOff) {
            // Get corresponding section numbers from spherical partition to classify its state
            secNum = positionOrientationPart->getSectionNumber(relativePosition, relativeOrientation, quatReference);
            discreteState  = maxNumberBoundStates + secNum;
        }
        // If none of the statements before modified discreteState, it returns the unbound state (0)
        return discreteState;
    };




    /* Auxiliary function used by sampleDiscreteState. Given two particles, use their positions and
     * orientations to determine if they are in one of the bound states. If not, return -1 (the coreMSM
     * approach can later assign the value of the previous state.) */
    template<int numBoundStates>
    int discreteTrajectory<numBoundStates>::getBoundState(vec3<double> relativePosition,
                                                          quaternion<double> relativeOrientation) {

        /* Check if it matches a bound states, if so return the corresponding state. Otherwise
         * return -1. */
        vec3<double> relPosCenter;
        quaternion<double> relQuatCenter;
        double angleDistance;
        for (int i = 0; i < boundStates.size(); i++) {
            relPosCenter = std::get<0>(boundStates[i]);
            relQuatCenter = std::get<1>(boundStates[i]);

            if ( (relPosCenter - relativePosition).norm() <= tolerancePosition) {
                angleDistance = msmrdtools::quaternionAngleDistance(relQuatCenter, relativeOrientation);
                if  ( angleDistance < toleranceOrientation) {
                    return i + 1;
                }
            }
        }
        return -1;
    };


    /*
     * Mostly only used when interacting using python and/or interacting with pybind.
     */

    /* From a given trajectory of the from (timestep, position, orientation), where repeated timesteps mean
     * differente particles at same tieme step, obtain a discrete trajectory using the patchyDimer discretization.
     * This is useful to load trajectories directly from python and discretize them.*/
    template<int numBoundStates>
    std::vector<double> discreteTrajectory<numBoundStates>::discretizeTrajectory(
            std::vector<std::vector<double>> trajectory) {
        int numParticles = 2; // Must be two to discretize trajectory (also it is a dimer)
        int timesteps = static_cast<int>(trajectory.size() / numParticles);

        // Set output trajectory
        std::vector<double> discreteTrajectory(timesteps);

        vec3<double> position1;
        vec3<double> position2;
        quaternion<double> orientation1;
        quaternion<double> orientation2;

        int prevDiscreteState = 0;
        int discreteState = 0;

        for (int i = 0; i < timesteps; i++) {
            auto part1Data = trajectory[numParticles*i];
            auto part2Data = trajectory[numParticles*i + 1];
            position1 = {part1Data[1], part1Data[2], part1Data[3]};
            position2 = {part2Data[1], part2Data[2], part2Data[3]};
            orientation1 = {part1Data[4], part1Data[5], part1Data[6], part1Data[7]};
            orientation2 = {part2Data[4], part2Data[5], part2Data[6], part2Data[7]};
            auto dummyParticle1 = particle(0, 0, position1, orientation1);
            auto dummyParticle2 = particle(0, 0, position2, orientation2);
            discreteState = this->sampleDiscreteState(dummyParticle1, dummyParticle2);
            // If sampleDiscreteState returned -1, return previous sample (CoreMSM approach).
            if (discreteState == -1) {
                discreteState = 1 * prevDiscreteState;
            }
            prevDiscreteState = 1*discreteState;

            discreteTrajectory[i] = discreteState;
        }
        return discreteTrajectory;
    }



    /* From a given trajectory H5 file of the from (timestep, position, orientation), where repeated timesteps mean
     * differente particles at same tieme step, obtain a discrete trajectory using the discreteTrajectory discretization.
     * This is the same as discretizeTrajectory, but loads the H5 file directly in c++ and later discretizes them.*/
    template<int numBoundStates>
    std::vector<double> discreteTrajectory<numBoundStates>::discretizeTrajectoryH5(std::string filename) {
        int numParticles = 2; // Must be two to discretize trajectory (also it is a dimer)
        vec3<double> position1;
        vec3<double> position2;
        quaternion<double> orientation1;
        quaternion<double> orientation2;

        int prevDiscreteState = 0;
        int discreteState = 0;

        // Read H5 file
        const H5std_string  FILE_NAME(filename);
        const H5std_string  DATASET_NAME("msmrd_data");
        H5File file(FILE_NAME, H5F_ACC_RDONLY);
        DataSet dataset = file.openDataSet(DATASET_NAME);

        // Get dimensions of dataset
        DataSpace dataspace = dataset.getSpace();
        hsize_t dims[2];
        int rank = dataspace.getSimpleExtentDims(dims);

        // Define output of read data
        int NX = dims[0];
        int NY = dims[1];
        double *trajectory = new double[NX*NY];

        // Set output trajectory
        int timesteps = static_cast<int>(dims[0] / numParticles);
        std::vector<double> discreteTrajectory(timesteps);

        // Define the memory space to read dataset.
        DataSpace mspace(rank, dims);

        dataset.read(trajectory, PredType::NATIVE_DOUBLE, mspace, dataspace);

        file.close();

        for (int i = 0; i < timesteps; i++) {
            /* 2D array indexes (row, col), with row = numParticles*i and col the column between 0 and 7. Its
             * 1D version index should be row*NY + col*/
            int jj = (numParticles * i)*NY;
            int kk = (numParticles*i + 1)*NY;
            position1 = {trajectory[jj+1], trajectory[jj+2], trajectory[jj+3]};
            position2 = {trajectory[kk+1], trajectory[kk+2], trajectory[kk+3]};
            orientation1 = {trajectory[jj+4], trajectory[jj+5], trajectory[jj+6], trajectory[jj+7]};
            orientation2 = {trajectory[kk+4], trajectory[kk+5], trajectory[kk+6], trajectory[kk+7]};
            auto dummyParticle1 = particle(0, 0, position1, orientation1);
            auto dummyParticle2 = particle(0, 0, position2, orientation2);
            discreteState = this->sampleDiscreteState(dummyParticle1, dummyParticle2);
            // If sampleDiscreteState returned -1, return previous value (CoreMSM approach).
            if (discreteState == -1) {
                discreteState = 1 * prevDiscreteState;
            }
            prevDiscreteState = 1*discreteState;

            discreteTrajectory[i] = discreteState;
        }
        delete [] trajectory;
        return discreteTrajectory;
    }


    /*
     * Setter functions to modify parameters
     */

    // Sets radial bounds of discretization transition region
    template<int numBoundStates>
    void discreteTrajectory<numBoundStates>::setRadialBounds(double rlower, double rupper) {
        rLowerBound = rlower;
        rUpperBound = rupper;
    };

    template<int numBoundStates>
    void discreteTrajectory<numBoundStates>::setTolerances(double positionTolerance, double orientationTolerance){
        tolerancePosition = positionTolerance;
        toleranceOrientation = orientationTolerance;
    };


}