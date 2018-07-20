//
// Created by dibakma on 22.06.18.
//

#pragma once
#include "quaternion.hpp"
#include "vec3.hpp"
//#include <pybind11/pybind11.h>
//#include <pybind11/numpy.h>
//#include <pybind11/stl.h>


/**
 * Declaration of the base class for particles
 */

class particle {
public:
    int pid;
    int type;
    double D;
    double Drot;
    vec3<double> position;
    quaternion<double> orientation;
    /**
     * Constructors
     * @param pid ID of the particle
     * @param type particle type, corresponds to msmid
     * @param D Diffusion constant
     * @param Drot rotational diffusion constant
     * @param position initial position of the particle
     * @param u normalized unit vector representing the initial orientation of the particle
     */
    particle(int pid, int type, double D, double Drot, vec3<double> position, quaternion<double> orientation):
            pid(pid), type(type), D(D), Drot(Drot), orientation(orientation), position(position){};

    particle(int pid, int type, double D, double Drot, std::vector<double> &position, std::vector<double> &u): pid(pid), type(type), D(D), Drot(Drot), position(position) {
        quaternion<double> q(0, u[0], u[1], u[2]);
        orientation = q;
    }

    /** Get properties functions for pybinding **/
    int getID() { return  pid; }
    int getType() { return  type; }
    int getD() { return  D; }
    int getDrot() { return  Drot; }
//    vec3<double> getPosition() { return  position; }
//    quaternion<double> getOrientation() { return  orientation; }
};