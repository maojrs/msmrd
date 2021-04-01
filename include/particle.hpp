//
// Created by dibakma on 22.06.18.
//

#pragma once
#include <map>
#include "quaternion.hpp"
#include "tools.hpp"
#include "vec3.hpp"

namespace msmrd {
    /**
     * Declaration of the class for particles. One line implementations in headers; otherwise in particle.cpp.
     */

    class particle {
    public:
        // Main variables
        int pid = 0;
        bool active = true;
        int type = 0;
        double D;
        double Drot;
        vec3<double> position;
        vec3<double> orientvector;
        quaternion<double> orientation;
        vec3<double> nextPosition;
        vec3<double> nextOrientvector;
        quaternion<double> nextOrientation;
        /* Next variables exclusive for Markovian switching behavior
         * (have default values if not initialized) */
        int state = 0;
        int nextState;
        double lagtime = 0;
        double timeCounter = 0;
        bool propagateTMSM = true;
        bool activeMSM = true;
        int boundTo = -1;
        int boundState = -1;
        /* Next variables exclusive for multiparticle or multiparticle MSM/RD behavior
         * (have default values if not initialized) */
        std::vector<int> boundList;
        std::vector<int> boundStates;
        int compoundIndex = -1;
        std::vector<int> activePatchList;
        /**
         * @param pid ID of the particle (can be the index in a particle list if it is part of particle list.)
         * @param active determines if particle currently active
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
         *
         * Next variables only for Markovian switching behavior
         * @param state particle current unbound state (begins in zero). If -1, then particle is in bound state.
         * @param nextState particle next state given and changed by the msm/ctmsm. Only useful in certain contexts.
         * @param lagtime saves the current lagtime from the MSM.
         * @param timeCounter counts time in between CTMSM iterations. It keeps track of time elapsed within one
         * timestep of the accumulated lagtimes from transitions of the CTMSM. Also used to store time left until
         * next reaction (activation) in the MAPK implementation.
         * @param propagateTMSM determines if MSM/CTMSM should be propagated on a given timestep
         * tcount and propagateTMSM are used to synchronize TMSM with diffusion/rotation timestepping
         * @param activeMSM determines if MSM behavior is active or dormant. This applied only to the MSM for
         * unbound particles. If in bound state in MSM/RD, this activeMSM should be set to false.
         * @param boundTo determines if particle is bound and to which particle in the particle list. If -1 the
         * particle is not bound, otherwise the value of boundTo is the index of the other particle in the particle
         * list.
         * @param boundState determines the state of the particle if it is bound to another particle. If it is not
         * bound, boundState = -1 and the state is given by the 'state' variable.
         *
         * Next variables are exclusive and only active for multiparticle MSM/RD functionality
         * @param boundList alternative to boundTo for multiparticle bindings. If empty, particle is not bound.
         * Otherwise the indexes stored in the vector correspond to the indices of the particles to which it is
         * bound in the particle list.
         * @param boundStates alternative to boundState for multiparticle bindings. If empty, particle is not bound.
         * Otherwise, the values correspond to the bound states of each of the corresponding particles in the boundList.
         * The indexes should be such that it is the same index as in the boundList.
         * @param compoundIndex if the particle is part of a particle compund; it corresponds to the index
         * in the particleCompounds vector in the multiparticle MSM/RD integrator. If it is -1, it means either
         * the particle does not belong to any particle compound, or that this functionality is not in use.
         * @param activePatchList If using patchy particles, this is a list of int values. The size
         * of the list correspond to the number of patches, a -1 value in the list indicates the corresponding
         * patch is unbound and active, a positive value corresponds to the index of the particle to which it
         * is bound, and it indicates the patch is only allowed to interact with the particle with which it is bound.
         * Only used for some multiparticle simulations.
         */

        // Constructors: receive input from vec3/quaternion or std::vector and numpy arrays (through pybind)
        particle(double D, double Drot, vec3<double> position, quaternion<double> orientation);

        particle(double D, double Drot, std::vector<double> &position, std::vector<double> &orientation);

        /* Constructors of particle with Markovian switching activated. Receive input from vec3/quaternion or
         * std::vector and numpy arrays (through pybind) */
        particle(int type, int state, double D, double Drot, vec3<double> position,
                 quaternion<double> orientation);

        particle(int type, int state, double D, double Drot, std::vector<double> &position,
                 std::vector<double> &orientation);

        /* The functions below of the particle class are used by c++ and python,
         * some of them are only useful to generate the python bindings. Functions
         * longer than one line are implemented in particle.cpp*/

        /* Getter functions */

        int getID() const { return pid; }

        double getD() const { return D; }

        double getDrot() const { return Drot; }

        int getType() const { return type; }

        // For Markovian switching behavior

        int getState() const { return state; }

        double getLagtime() const { return lagtime; }

        int getBoundTo() const { return boundTo; }

        int getBoundState() const { return boundState; }

        std::vector<int> getBoundList() const { return boundList; }

        std::vector<int> getBoundStates() const { return boundStates; }

        int getCompoundIndex() const { return compoundIndex; }



        // Setter functions

        void setID(int newID) { pid = newID; }

        void setD(double Dnew) { D = Dnew; }

        void setDrot(double Drotnew) { Drot = Drotnew; }

        void setDs(double Dnew, double Drotnew) { D = Dnew, Drot = Drotnew; }

        void setType(int newtype) { type = newtype; }

        void setPosition(vec3<double> newposition) { position = newposition; }

        void setNextPosition(vec3<double> nextposition) { nextPosition = nextposition; }

        void setOrientVector(vec3<double> neworientvector) { orientvector = neworientvector; }

        void setNextOrientVector(vec3<double> nextorientvector) { nextOrientvector = nextorientvector; }

        void setOrientation(quaternion<double> neworientation) { orientation = neworientation; }

        void setNextOrientation(quaternion<double> nextorientation) { nextOrientation = nextorientation; }

        // Exclusive for pybindings

        void setPositionPyBind(std::vector<double> newposition) { position = newposition; }

        void setOrientVectorPyBind(std::vector<double> neworientvector) { orientvector = neworientvector; }

        void setOrientationPyBind(std::vector<double> neworientation);

        // For Markovian switching behavior

        void setState(int newstate);

        void setNextState(int nextstate) { nextState = nextstate; }

        void setLagtime(double newlagtime) { lagtime = newlagtime; }

        void setBoundTo(int particleIndex) { boundTo = particleIndex; }

        void setBoundState(int newBoundState);

        void setMSMoff() { activeMSM = false; }

        void setMSMon() { activeMSM = true; }

        /* Additional functions. */
        void updatePosition();

        void updateOrientation();

        void activate() { active = true; }

        void deactivate() { active = false; }

        bool isActive() { return active; }

        // For Markovian switching behavior

        void updateState() { state = 1 * nextState; };

        bool isMSMactive() { return activeMSM; }

        void activateMSM() { activeMSM = true; }

        void deactivateMSM() { activeMSM = false; }

        void deactivateResetMSM();

        void setActivePatchList(int numPatches);

    };


}




