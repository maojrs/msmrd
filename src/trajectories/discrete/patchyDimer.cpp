//
// Created by maojrs on 1/22/19.
//
#include "trajectories/discrete/patchyDimer.hpp"

namespace msmrd {

    /*
     * Chooses how to discretize the full trajectory (trajectoryData) into a discretized trajectory
     * to be analyzed and extracted into a Markov state model. This specific discretization will follow
     * the core MSM approach
     */
    patchyDimer::patchyDimer(unsigned long Nparticles, int bufferSize) :
    trajectoryPositionOrientation(Nparticles, bufferSize) {
        int numSphericalSectionsPos = 7; //7; //7;
        int numRadialSectionsQuat = 5; // 3; //5;
        int numSphericalSectionsQuat =7; // 6; //7;
        discreteTrajectoryData.reserve(bufferSize);
        positionOrientationPart = std::make_unique<positionOrientationPartition>(rUpperBound, numSphericalSectionsPos,
                                                                                 numRadialSectionsQuat,
                                                                                 numSphericalSectionsQuat);
        setMetastableRegions();
    };

    patchyDimer::patchyDimer(unsigned long Nparticles, int bufferSize, double rLowerBound, double rUpperBound) :
        trajectoryPositionOrientation(Nparticles, bufferSize), rLowerBound(rLowerBound), rUpperBound(rUpperBound) {
        int numSphericalSectionsPos = 7; //7; //7;
        int numRadialSectionsQuat = 5; // 3; //5;
        int numSphericalSectionsQuat =7; // 6; //7;
        discreteTrajectoryData.reserve(bufferSize);
        positionOrientationPart = std::make_unique<positionOrientationPartition>(rUpperBound, numSphericalSectionsPos,
                                                                                 numRadialSectionsQuat,
                                                                                 numSphericalSectionsQuat);
        setMetastableRegions();
    };

    /* Samples discrete trajectory directly from the particle list. This is done for every timestep and saved into
     * the class data directly, so this function is better suited to obtain the discrete trajectory at the same
     * time as it is being computed. It is assumed only the first two particles are relevant. */
    void patchyDimer::sampleDiscreteTrajectory(double time, std::vector<particle> &particleList) {
        // Initialize sample with value zero
        int sample = sampleDiscreteState(particleList[0], particleList[1]);

        // Save previous value and push into trajectory
        prevsample = 1*sample;
        discreteTrajectoryData.push_back(std::vector<int>{sample});
    };


    /* Main function to sample the discrete state of two particles. It returns the corresponding
     * bound state, transition state or unbound state (0). In the bound region (r< rLowerBound), it uses
     * the core MSM approach to assign a value (previous state if current state is not a bound or
     * transition state). */
    int patchyDimer::sampleDiscreteState(particle part1, particle part2) {
        // Initialize sample with value zero
        int discreteState = 0;

        /* Calculate relative position taking into account periodic boundary measured
         * from i to j (gets you from i to j). */
        vec3<double> relativePosition = calculateRelativePosition(part1.position, part2.position);

        // Rotate relative position to match the reference orientation of particle 1. (VERY IMPORTANT)
        relativePosition = msmrdtools::rotateVec(relativePosition, part1.orientation.conj());
        quaternion<double> quatReference = {1,0,0,0}; // we can then define reference quaternion as identity.

        // Calculate relative orientation (w/respect to particle 1)
        quaternion<double> relativeOrientation;
        relativeOrientation = part1.orientation.conj() * part2.orientation;

        // Extract current state, save into sample and return sample
        int secNum;
        if (relativePosition.norm() < rLowerBound) {
            discreteState = getBoundState(relativePosition, relativeOrientation);
            // If sample doesn't correspond to any bound state in r<rLowerBound, assign previous state (CoreMSM)
            if (discreteState == -1) {
                    discreteState = prevsample;
            }
        } else if (relativePosition.norm() < positionOrientationPart->relativeDistanceCutOff) {
            // Get corresponding section numbers from spherical partition to classify its state
            secNum = positionOrientationPart->getSectionNumber(relativePosition, relativeOrientation, quatReference);
            discreteState  = maxNumberBoundStates + secNum;
        }
        return discreteState;
        //discreteTrajectoryData.push_back(sample);
    };


    /* Given two particles, use their positions and orientations to determine if they are in one of
     * the 8 bound states (1 to 8). If not, return -1 (to later assign the value of the previous state) */
    int patchyDimer::getBoundState(vec3<double> relativePosition, quaternion<double> relativeOrientation) {

        /* Check if it matches a bound states, if so return the corresponding state. Otherwise
         * return -1. */
        vec3<double> relPosCenter;
        quaternion<double> relQuatCenter;
        double angleDistance;
        for (int i = 0; i < 8; i++) {
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


    /* Sets metastable regions (bound states) of patchy dimer (two equal patcy particles, each with two patches
     * an ang angleDiff away). The centers of the metastable regions are given by a tuple of relative position
     * and relative orientation. The size of the regions are determined by tolerancePosition and
     * toleranceOrientation*/
    void patchyDimer::setMetastableRegions() {
        double angleDiff = 3 * M_PI / 5; // angle difference to form a pentamer
        /* Define relative position vectors from particle 1 at the origin. These two patches
         * point in the same direction as the two patches in the dimer. */
        vec3<double> relPos1 = {std::cos(angleDiff / 2.0), std::sin(angleDiff / 2.0), 0};
        vec3<double> relPos2 = {std::cos(angleDiff / 2.0), std::sin(-angleDiff / 2.0), 0};
        vec3<double> relPos1orthogonal = {-1.0 * std::sin(angleDiff / 2.0), std::cos(angleDiff / 2.0), 0.0};
        vec3<double> relPos2orthogonal = {std::sin(angleDiff / 2.0), std::cos(angleDiff / 2.0), 0.0};
        /* Relative rotations (from particle 1) of particle 2 that yield the 8 bound states
         * in the axis-angle representation. (One needs to make drawing to understand)*/
        std::array<vec3<double>, 8> rotations;
        rotations[0] = M_PI * relPos1orthogonal; //ok
        rotations[1] = {0.0, 0.0, -2 * M_PI / 5.0}; //ok
        rotations[2] = {0.0, 0.0, M_PI}; //ok
        rotations[3] = {0.0, M_PI, 0.0}; //ok
        // --first 4 rotations correspond to binding on top patch of particle 1, next 4 rotations to bottom patch
        rotations[4] = M_PI * relPos2orthogonal; //ok
        rotations[5] = {0.0, 0.0, 2 * M_PI / 5.0}; //ok
        rotations[6] = {0.0, 0.0, M_PI}; //ok
        rotations[7] = {0.0, M_PI, 0.0}; //ok
        /*Convert rotations in the axis angle representation to quaternions */
        std::array<quaternion<double>, 8> quatRotations;
        for (int i = 0; i < 8; i++) {
            quatRotations[i] = msmrdtools::axisangle2quaternion(rotations[i]);
        }
        // Fill bound states with corresponding combinations of relative position vectors and quaternion orientations.
        boundStates[0] = std::make_tuple(relPos1, quatRotations[0]);
        boundStates[1] = std::make_tuple(relPos1, quatRotations[1]);
        boundStates[2] = std::make_tuple(relPos1, quatRotations[2]);
        boundStates[3] = std::make_tuple(relPos1, quatRotations[3]);
        boundStates[4] = std::make_tuple(relPos2, quatRotations[4]);
        boundStates[5] = std::make_tuple(relPos2, quatRotations[5]);
        boundStates[6] = std::make_tuple(relPos2, quatRotations[6]);
        boundStates[7] = std::make_tuple(relPos2, quatRotations[7]);
    }


    /*
     * Mostly only used when interacting using python and/or interacting with pybind.
     */

    /* From a given trajectory of the from (timestep, position, orientation), where repeated timesteps mean
     * differente particles at same tieme step, obtain a discrete trajectory using the patchyDimer discretization.
     * This is useful to load trajectories directly from python and discretize them.*/
    std::vector<double> patchyDimer::discretizeTrajectory(std::vector<std::vector<double>> trajectory) {
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
            discreteState = getState(dummyParticle1, dummyParticle2);
            // If getState returned -1, return previous (CoreMSM approach).
            if (discreteState == -1) {
                discreteState = 1 * prevDiscreteState;
            }
            prevDiscreteState = 1*discreteState;

            discreteTrajectory[i] = discreteState;
        }
        return discreteTrajectory;
    }



    /* From a given trajectory H5 file of the from (timestep, position, orientation), where repeated timesteps mean
     * differente particles at same tieme step, obtain a discrete trajectory using the patchyDimer discretization.
     * This is the same as discretizeTrajectory, but loads the H5 file directly in c++ and later discretizes them.*/
    std::vector<double> patchyDimer::discretizeTrajectoryH5(std::string filename) {
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
        //double trajectory[320000][8];
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
            discreteState = getState(dummyParticle1, dummyParticle2);
            // If getState returned -1, return previous (CoreMSM approach).
            if (discreteState == -1) {
                discreteState = 1 * prevDiscreteState;
            }
            prevDiscreteState = 1*discreteState;

            discreteTrajectory[i] = discreteState;
        }
        delete [] trajectory;
        return discreteTrajectory;
    }


    /* Similar to getBoundState, but it returns the corresponding bound state, transition state or
     * unbound state (0). If it returns -1, then the trajectory is below the rLowerBound, but it
     * does not belong to any bound state; therefore, the previous state sould be assigned if computing a
     * discrete trajectory. It is a mixture between sampleDiscreteState and getBoundState. It will
     * specially useful for PyBind when calculating benchmarks in python interface and when using
     * discretizeTrajectory to obtain discrete trajectories directly from python arrays.  */
    int patchyDimer::getState(particle part1, particle part2) {
        // Calculate relative distance taking into account periodic boundary.
        vec3<double> relativePosition = calculateRelativePosition(part1.position, part2.position);

        // Rotate relative position to match the reference orientation of particle 1.
        relativePosition = msmrdtools::rotateVec(relativePosition, part1.orientation.conj());
        quaternion<double> quatReference = {1,0,0,0}; // we can then define reference quaternion as identity.

        // Calculate relative orientation (w/respect to particle 1)
        quaternion<double> relativeOrientation = part1.orientation.conj() * part2.orientation;

        // Check if it matches a bound state, if so return the corresponding state. Otherwise return -1.
        int secNum;
        int state;
        vec3<double> relPosCenter;
        quaternion<double> relQuatCenter;
        double angleDistance;
        // Returns bound state, -1 or transitionState if below rLowerBound region but not in any bound state
        if (relativePosition.norm() < rLowerBound) {
            for (int i = 0; i < 8; i++) {
                relPosCenter = std::get<0>(boundStates[i]);
                relQuatCenter = std::get<1>(boundStates[i]);
                // If in bound state, returns bound state
                if ((relPosCenter - relativePosition).norm() <= tolerancePosition) {
                    angleDistance = msmrdtools::quaternionAngleDistance(relQuatCenter, relativeOrientation);
                    if (angleDistance < toleranceOrientation) {
                        return i + 1; // returns bound state
                    }
                }
            }
            return -1; // returns -1 when not in bound state but in r<rLowerBound, so CoreMSM is later applied.
        }
        // Returns transition state
        else if (relativePosition.norm() < rUpperBound) {
            // Get corresponding section numbers from spherical partition to classify its state
            secNum = positionOrientationPart->getSectionNumber(relativePosition, relativeOrientation, quatReference);
            return maxNumberBoundStates + secNum;
        }
        // Returns unbound state
        else {
            return 0;
        }
    };

}