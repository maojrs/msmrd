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
    std::string bodytype;
    vec3<double> position;
    vec3<double> orientvector;
    quaternion<double> orientation;
    /**
     * @param pid ID of the particle
     * @param active determines if particle currently active
     * @param D diffusion constant
     * @param Drot rotational diffusion constant
     * @param bodytype determines rotation integrator behavior, can be either point, rod or rigidsolid
     * (determined by orientational degrees of freedom, points have no orientation, rods need only one vector and
     * rigidsolid requires a complete quaternion).
     * @param position initial position of the particle
     * @param orientvector orientation vector (only to be used by for rod-like particles)
     * @param orientation normalized quaternion representing the initial orientation of the particle
     */

    // Constructors: receive input from vec3/quaternion or std::vector and numpy arrays (through pybind)
    particle(double D, double Drot, std::string bodytype, vec3<double> position, quaternion<double> orientation)
            : D(D), Drot(Drot), bodytype(bodytype), position(position), orientation(orientation){
        orientvector = vec3<double>(0.,0.,1.);
    };

    particle(double D, double Drot, std::string bodytype, std::vector<double> &position, std::vector<double> &orientation)
            : D(D), Drot(Drot), bodytype(bodytype), position(position), orientation(orientation) {
        orientvector = vec3<double>(0.,0.,1.);
    };

    /**
     * Get and set functions. Some used by c++ and python,
     * some only to be used by pyhon with python bindings.
     **/
    int getID() { return  pid; }
    double getD() { return  D; }
    double getDrot() { return  Drot; }
    std::string getBodyType() { return  bodytype; }
    void setD(double Dnew) { D = Dnew; }
    void setDrot(double Drotnew) { Drot = Drotnew; }
    void setBodyType(std::string bodytypenew) { bodytype = bodytypenew; }
    void setPosition(vec3<double> newposition) { position = newposition; }
    void setPositionPyBind(std::vector<double> newposition) { position = newposition; }
    void setOrientVector(vec3<double> neworientvector) { orientvector = neworientvector; }
    void setOrientVectorPyBind(std::vector<double> neworientvector) { orientvector = neworientvector; }
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
    particleMS(int type, int state, double D, double Drot, std::string bodytype, vec3<double> position, quaternion<double> orientation)
            : type(type), state(state), particle(D, Drot, bodytype, position, orientation){};

    particleMS(int type, int state, double D, double Drot, std::string bodytype, std::vector<double> &position, std::vector<double> &orientation)
            : type(type), state(state), particle(D, Drot, bodytype, position, orientation) {};


    // Additional get and set functions for particleMS
    int getType() { return  type; }
    int getState() { return  state; }
    double getLagtime() { return  lagtime; }
    void setState(int newstate) { state = newstate; }
    void setType(int newtype) { type = newtype; }
    void setLagtime(double newlagtime) { lagtime = newlagtime; }
};
