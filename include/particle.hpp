//
// Created by dibakma on 22.06.18.
//

#pragma once
#include "quaternion.hpp"
#include "vec3.hpp"

/**
 * Declaration of the class for particles
 */
class particle {
protected:
    int pid = 0;
    bool active = true;
    vec3<double> nextPosition;
    vec3<double> nextOrientvector;
    quaternion<double> nextOrientation;
    /**
     * @param pid ID of the particle
     * @param active determines if particle currently active
     * @param nextPosition saves next position for integrator to update
     * @param nextOrientvector saves next orientation vector for integrator to update
     * @param nextOrientation saves next orientation quaternions for integrator to update
     */
public:
    double D;
    double Drot;
    std::string bodytype;
    vec3<double> position;
    vec3<double> orientvector;
    quaternion<double> orientation;
    /**
     * @param D diffusion constant
     * @param Drot rotational diffusion constant
     * @param bodytype determines rotation integrator behavior, can be either point, rod or rigidsolid
     * (determined by orientational degrees of freedom, points have no orientation, rods need only one vector and
     * rigidsolid requires a complete quaternion).
     * @param position position vector of the particle
     * @param orientvector orientation vector (only to be used by for rod-like particles)
     * @param orientation normalized quaternion representing the initial orientation of the particle
     */

    // Constructors: receive input from vec3/quaternion or std::vector and numpy arrays (through pybind)
    particle(double D, double Drot, std::string bodytype, vec3<double> position, quaternion<double> orientation)
            : D(D), Drot(Drot), bodytype(bodytype), position(position), orientation(orientation){
        orientvector = vec3<double>(0.,0.,1.);
        nextPosition = 1.0*position;
        nextOrientation = 1.0*orientation;
        nextOrientvector = 1.0*orientvector;
    };

    particle(double D, double Drot, std::string bodytype, std::vector<double> &position, std::vector<double> &orientation)
            : D(D), Drot(Drot), bodytype(bodytype), position(position), orientation(orientation) {
        orientvector = vec3<double>(0.,0.,1.);
        std::vector<double> nextPosition(position);
        std::vector<double> nextOrientation(orientation);
        nextOrientvector = 1.0*orientvector;
    };

    /* Additional functions and getters and setters. Some used by c++ and python,
     * some only to be used by pyhon with python bindings. */
    void updatePosition() { position = 1*nextPosition; };
    void updateOrientation() {
        orientvector = 1*nextOrientvector;
        orientation = 1*nextOrientation;
    };
    int getID() const { return  pid; }
    double getD() const { return  D; }
    double getDrot() const { return  Drot; }
    std::string getBodyType() const { return  bodytype; }

    void setD(double Dnew) { D = Dnew; }
    void setDrot(double Drotnew) { Drot = Drotnew; }
    void setBodyType(std::string bodytypenew) { bodytype = bodytypenew; }

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
};

/**
 * Declaration of particles subclass with Markovian Switch (MS)
 */
class particleMS: public particle {
protected:
    int nextState;
public:
    int type;
    int state;
    double lagtime = 0;
    double tcount = 0;
    double propagateTMSM = true ;
    /**
     * @param type particle type, corresponds to msmid
     * @param state particle current state
     * @param nextState particle next state given and changed by the msm/ctmsm
     * @param lagtime saves the current lagtime from the MSM
     * @param tcount counts time in between MSM/CTMSM iterations
     * @param propagateTMSM determines if MSM/CTMSM should be propagated on a given timestep
     * tcount and propagateTMSM are used to synchronize TMSM with diffusion/rotation timestepping
     */

    // Constructors: receive input from vec3/quaternion or std::vector and numpy arrays (through pybind)
    particleMS(int type, int state, double D, double Drot, std::string bodytype, vec3<double> position, quaternion<double> orientation)
            : type(type), state(state), particle(D, Drot, bodytype, position, orientation){};

    particleMS(int type, int state, double D, double Drot, std::string bodytype, std::vector<double> &position, std::vector<double> &orientation)
            : type(type), state(state), particle(D, Drot, bodytype, position, orientation) {};

    // Additional functions and getters and setters for particleMS
    void updateState() { state = 1*nextState; };
    int getType() const { return  type; }
    int getState() const { return  state; }
    double getLagtime() const { return  lagtime; }
    void setState(int newstate) { state = newstate; }
    void setNextState(int newstate) { nextState = newstate; }
    void setType(int newtype) { type = newtype; }
    void setLagtime(double newlagtime) { lagtime = newlagtime; }
};






