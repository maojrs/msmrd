//
// Created by dibakma on 22.06.18.
//

#pragma once
#include "quaternion.h"
#include "vec3.h"

/**
 * Declaration of the base class for particles
 * @tparam scalar underlying data type of the position and orientation
 */

template<typename scalar>
class particle {
public:
    int pid;
    int type;
    quaternion<scalar> orientation;
    vec3<scalar> position;
    scalar D;
    scalar Drot;

    particle(int pid, int type, scalar D, scalar Drot, quaternion<scalar> orientation, vec3<scalar> position):
            pid(pid), type(type), D(D), Drot(Drot), orientation(orientation), position(position){};

    /**
     * Specific constructor for double type
     * @param pid ID of the particle
     * @param type particle type
     * @param D Diffusion constant
     * @param Drot rotational diffusion constant
     * @param u normalized unit vector representing the initial orientation of the particle
     * @param position initial position of the particle
     */
    particle(int pid, int type, scalar D, scalar Drot, std::vector<double> &u, std::vector<double> &position): pid(pid), type(type), D(D), Drot(Drot), position(position) {
        quaternion<double> q(0, u[0], u[1], u[2]);
        orientation = q;
    }
};