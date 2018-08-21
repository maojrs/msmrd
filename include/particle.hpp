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
public:
    double D;
    double Drot;
    vec3<double> position;
    quaternion<double> orientation;
    /**
     * @param pid ID of the particle
     * @param active determines if particle currently active
     * @param D diffusion constant
     * @param Drot rotational diffusion constant
     * @param position initial position of the particle
     * @param orientation normalized quaternion representing the initial orientation of the particle
     */

    // Constructors: receive input from vec3/quaternion or std::vector and numpy arrays (through pybind)
    particle(double D, double Drot, vec3<double> position, quaternion<double> orientation)
            : D(D), Drot(Drot), position(position), orientation(orientation){};

    particle(double D, double Drot, std::vector<double> &position, std::vector<double> &orientation)
            : D(D), Drot(Drot), position(position), orientation(orientation) {};

    /**
     * Get and set functions. Some used by c++ and python,
     * some only to be used by pyhon with python bindings.
     **/
    int getID() { return  pid; }
    double getD() { return  D; }
    double getDrot() { return  Drot; }
    void setD(double Dnew) { D = Dnew; }
    void setDrot(double Drotnew) { Drot = Drotnew; }
    void setPosition(vec3<double> newposition) { position = newposition; }
    void setPositionPyBind(std::vector<double> newposition) { position = newposition; }
    void setOrientation(quaternion<double> neworientation) { orientation = neworientation; }
    void setOrientationPyBind(std::vector<double> neworientation) {
        quaternion<double> quat(neworientation);
        orientation = quat;
    }
};


/**
 * Declaration of particles subclass with Markovian Switch (MS)
 */
class particleMS: public particle {
public:
    int type;
    int state;
    double lagtime = 0;
    double tcount = 0;
    double propagateTMSM = true ;
    /**
     * @param type particle type, corresponds to msmid
     * @param state particle state given and changed by the msm/ctmsm
     * @param lagtime saves the current lagtime from the MSM
     * @param tcount counts time in between MSM/CTMSM iterations
     * @param propagateTMSM determines if MSM/CTMSM should be propagated on a given timestep
     * tcount and propagateTMSM are used to synchronize TMSM with diffusion/rotation timestepping
     */

    // Constructors: receive input from vec3/quaternion or std::vector and numpy arrays (through pybind)
    particleMS(int type, int state, double D, double Drot, vec3<double> position, quaternion<double> orientation)
            : type(type), state(state), particle(D, Drot, position, orientation){};

    particleMS(int type, int state, double D, double Drot, std::vector<double> &position, std::vector<double> &orientation)
            : type(type), state(state), particle(D, Drot, position, orientation) {};


    // Additional get and set functions for particleMS
    int getType() { return  type; }
    int getState() { return  state; }
    double getLagtime() { return  lagtime; }
    void setState(int newstate) { state = newstate; }
    void setType(int newtype) { type = newtype; }
    void setLagtime(double newlagtime) { lagtime = newlagtime; }
};
