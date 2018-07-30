//
// Created by dibakma on 22.06.18.
//

#pragma once
#include "quaternion.hpp"
#include "vec3.hpp"

/**
 * Declaration of the base class for particles
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
     * Constructors
     * @param pid ID of the particle
     * @param D Diffusion constant
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
    int getD() { return  D; }
    int getDrot() { return  Drot; }
    void setPosition(vec3<double> newposition) { position = newposition; }
    void setPositionPybind(std::vector<double> newposition) { position = newposition; }
    void setOrientation(quaternion<double> neworientation) { orientation = neworientation; }
    void setOrientationPybind(std::vector<double> neworientation) {
        quaternion<double> quat(neworientation);
        orientation = quat;
    }
};

class particleMS: public particle {
public:
    int type; //particle type, corresponds to msmid
    int state; //particle state given and changed by the msm

    // Constructors: receive input from vec3/quaternion or std::vector and numpy arrays (through pybind)
    particleMS(int type, int state, double D, double Drot, vec3<double> position, quaternion<double> orientation)
            : type(type), state(state), particle(D, Drot, position, orientation){};

    particleMS(int type, int state, double D, double Drot, std::vector<double> &position, std::vector<double> &orientation)
            : type(type), state(state), particle(D, Drot, position, orientation) {};

    int getType() { return  type; }
    int getState() { return  state; }
    void setState(int newstate) { state = newstate; }
    void setType(int newtype) { type = newtype; }


};
