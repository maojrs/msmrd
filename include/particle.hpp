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
        int type = 0;
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
         * @param type particle type (defaults to zero), for Markovian switching should be
         * considered the same as the msmid.
         * @param D diffusion constant
         * @param Drot rotational diffusion constant
         * @param position position vector of the particle
         * @param orientvector orientation vector (only to be used by for rod-like particles)
         * @param orientation normalized quaternion representing the initial orientation of the particle
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
            nextPosition = vec3<double>(position);
            nextOrientation = quaternion<double>(orientation);
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

        void setDs(double Dnew, double Drotnew) { D = Dnew, Drot = Drotnew; }

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
        int boundTo = -1;
        int boundState = -1;
        /**
         * @param state particle current unbound state. If -1, then particle is in bound state.
         * @param nextState particle next state given and changed by the msm/ctmsm
         * @param lagtime saves the current lagtime from the MSM
         * @param tcount counts time in between MSM/CTMSM iterations
         * @param propagateTMSM determines if MSM/CTMSM should be propagated on a given timestep
         * tcount and propagateTMSM are used to synchronize TMSM with diffusion/rotation timestepping
         * @param activeMSM determines if MSM behavior is active or dormant. This applied only to the MSM for
         * unbound particles. If in bound state, this activeMSM should be set to false.
         * @param boundTo determines if particle is bound and to which particle in the particle list. If -1 the
         * particle is not bound, otherwise the value of boundTo is the index of the other particle in the particle
         * list.
         * @param boundState determines the state of the particle if it is bound to another particle. If it is not
         * bound, boundState = -1 and the state is given by the 'state' variable.
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

        bool isMSMactive() {return activeMSM; }

        void activateMSM() {activeMSM = true; }

        void deactivateMSM() {activeMSM = false; }

        int getState() const { return state; }

        double getLagtime() const { return lagtime; }

        int getBoundTo() const { return boundTo; }

        int getBoundState() const { return boundState; }

        void setState(int newstate) {
            state = newstate;
            boundState = -1;
            boundTo = -1;
        }

        void setNextState(int nextstate) { nextState = nextstate; }

        void setLagtime(double newlagtime) { lagtime = newlagtime; }

        void setBoundTo(int particleIndex) { boundTo = particleIndex; }

        void setBoundState(int newBoundState) {
            state = -1;
            boundState = newBoundState;
            activeMSM = false;
        }

        void setNextType(int nexttype) { nextType = nexttype; }

        void setMSMoff() {activeMSM = false; }

        void setMSMon() {activeMSM = true; }
    };

}




