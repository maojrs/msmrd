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
public:
    int pid;
    int type;
    int state;
    double D;
    double Drot;
    vec3<double> position;
    quaternion<double> orientation;
    /**
     * Constructors
     * @param pid ID of the particle
     * @param type particle type, corresponds to msmid
     * @param state particle state given and changed by the msm
     * @param D Diffusion constant
     * @param Drot rotational diffusion constant
     * @param position initial position of the particle
     * @param orientation normalized quaternion representing the initial orientation of the particle
     */
    // Constructors: receive input from vec3/quaternion or std::vector and numpy arrays (through pybind)
    particle(int pid, int type, int state, double D, double Drot, vec3<double> position, quaternion<double> orientation):
            pid(pid), type(type), state(state), D(D), Drot(Drot), orientation(orientation), position(position){};

    particle(int pid, int type, int state, double D, double Drot, std::vector<double> &position, std::vector<double> &orientation)
            : pid(pid), type(type), state(state), D(D), Drot(Drot), position(position), orientation(orientation) {};

    /** Get properties functions for pybinding **/
    int getID() { return  pid; }
    int getType() { return  type; }
    int getState() { return  state; }
    int getD() { return  D; }
    int getDrot() { return  Drot; }
};