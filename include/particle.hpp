//
// Created by dibakma on 22.06.18.
//

#pragma once
#include "quaternion.hpp"
#include "tools.hpp"
#include "vec3.hpp"

namespace msmrd {
    /**
     * Declaration of the class for particles
     */
    class particle {
    protected:
        int pid = 0;
        bool active = true;
        /**
         * @param pid ID of the particle
         * @param active determines if particle currently active
         */
    public:
        int type;
        double D;
        double Drot;
        vec3<double> position;
        vec3<double> orientvector;
        quaternion<double> orientation;
        vec3<double> nextPosition;
        vec3<double> nextOrientvector;
        quaternion<double> nextOrientation;
        std::vector<quaternion<double>> symmetryQuaternions;
        /**
         * @param D diffusion constant
         * @param Drot rotational diffusion constant
         * @param position position vector of the particle
         * @param orientvector orientation vector (only to be used by for rod-like particles)
         * @param orientation normalized quaternion representing the initial orientation of the particle
         * @param type particle type (defaults to zero), for Markovian swithching should be
         * considered the same as the msmid.
         * @param nextPosition saves next position for integrator to update
         * @param nextOrientvector saves next orientation vector for integrator to update
         * @param nextOrientation saves next orientation quaternions for integrator to update
         * @param symmetryQuats list of quaternions that represent symmetry rotations for the particle.
         */

        // Constructors: receive input from vec3/quaternion or std::vector and numpy arrays (through pybind)
        particle(double D, double Drot, vec3<double> position, quaternion<double> orientation)
                : D(D), Drot(Drot), position(position), orientation(orientation) {
            type = 0;
            orientvector = vec3<double>(0., 0., 1.);
            orientvector = msmrdtools::rotateVec(orientvector, orientation);
            nextPosition = 1.0 * position;
            nextOrientation = 1.0 * orientation;
            nextOrientvector = 1.0 * orientvector;
            symmetryQuaternions.resize(0);
        };

        particle(double D, double Drot, std::vector<double> &position,
                 std::vector<double> &orientation)
                : D(D), Drot(Drot), position(position), orientation(orientation) {
            type = 0;
            orientvector = vec3<double>(0., 0., 1.);
            orientvector = msmrdtools::rotateVec(orientvector, orientation);
            std::vector<double> nextPosition(position);
            std::vector<double> nextOrientation(orientation);
            nextOrientvector = 1.0 * orientvector;
            symmetryQuaternions.resize(0);
        };

        /* Additional functions and getters and setters. Some used by c++ and python,
         * some only to be used by pyhon with python bindings. */
        void updatePosition() { position = 1.0 * nextPosition; };

        void updateOrientation() {
            orientvector = 1 * nextOrientvector;
            orientation = 1 * nextOrientation;
        };

        void activate() { active = true; }

        void deactivate() { active = false; }

        bool isActive() { return active; }

        int getID() const { return pid; }

        double getD() const { return D; }

        double getDrot() const { return Drot; }

        int getType() const { return type; }

        void setD(double Dnew) { D = Dnew; }

        void setDrot(double Drotnew) { Drot = Drotnew; }

        void setType(int newtype) { type = newtype; }

        void setPosition(vec3<double> newposition) { position = newposition; }

        void setNextPosition(vec3<double> nextposition) { nextPosition = nextposition; }

        void setPositionPyBind(std::vector<double> newposition) { position = newposition; }

        void setOrientVector(vec3<double> neworientvector) { orientvector = neworientvector; }

        void setNextOrientVector(vec3<double> nextorientvector) { nextOrientvector = nextorientvector; }

        void setOrientVectorPyBind(std::vector<double> neworientvector) { orientvector = neworientvector; }

        void setOrientation(quaternion<double> neworientation) { orientation = neworientation; }

        void setNextOrientation(quaternion<double> nextorientation) { nextOrientation = nextorientation; }

        void setOrientationPyBind(std::vector<double> neworientation) {
            quaternion<double> quat(neworientation);
            orientation = quat;
        }

        void setSymmetryQuaternions(int numQuats, std::vector<std::vector<double>> &symQuaternionsList) {
            symmetryQuaternions.resize(numQuats);
            for (int i = 0; i < numQuats; i ++){
                symmetryQuaternions[i] = symQuaternionsList[i];
            }
        }
    };

    /**
     * Declaration of particles subclass with Markovian Switch (MS)
     */
    class particleMS : public particle {
    protected:
        int nextType;
        int nextState;
    public:
        int state;
        double lagtime = 0;
        double tcount = 0;
        bool propagateTMSM = true;
        bool activeMSM = true;
        std::array<int, 2> indexesBoundParts{{-1, -1}};
        /**
         * @param state particle current state
         * @param nextState particle next state given and changed by the msm/ctmsm
         * @param lagtime saves the current lagtime from the MSM
         * @param tcount counts time in between MSM/CTMSM iterations
         * @param propagateTMSM determines if MSM/CTMSM should be propagated on a given timestep
         * tcount and propagateTMSM are used to synchronize TMSM with diffusion/rotation timestepping
         * @param activeMSM determines if MSM behavior is active or dormant. Under certain conditions (bound state),
         * maybe it is convenient to turn off the MSM behavior.
         * @param indexesBoundParts if particle is being used to model a bound state between two other particles,
         * this variable indicate the indexes of the two bound particles in the particleList being intergrated. If
         * not used, they take the default value -1.
         */

        // Constructors: receive input from vec3/quaternion or std::vector and numpy arrays (through pybind)
        particleMS(int type0, int state, double D, double Drot, vec3<double> position,
                   quaternion<double> orientation)
                : state(state), particle(D, Drot, position, orientation) {
            type = type0;
        };

        particleMS(int type0, int state, double D, double Drot, std::vector<double> &position,
                   std::vector<double> &orientation)
                : state(state), particle(D, Drot, position, orientation) {
            type = type0;
        };

        // Additional functions and getters and setters for particleMS
        void updateState() { state = 1 * nextState; };

        void updateType() { type = 1 * nextType; };

        int getState() const { return state; }

        double getLagtime() const { return lagtime; }

        void setState(int newstate) { state = newstate; }

        void setNextState(int nextstate) { nextState = nextstate; }

        void setNextType(int nexttype) { nextType = nexttype; }

        void setLagtime(double newlagtime) { lagtime = newlagtime; }

        void setMSMoff() {activeMSM = false; }

        void setMSMon() {activeMSM = true; }
    };

}




