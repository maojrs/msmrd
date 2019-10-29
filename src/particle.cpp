//
// Created by maojrs on 10/17/19.
//

#include "particle.hpp"

namespace msmrd {

    /**
     * Implementation of particle class
     */

    /* Constructors of normal particles: receive input from vec3/quaternion or std::vector and
     * numpy arrays (through pybind) */
    particle::particle(double D, double Drot, vec3<double> position, quaternion<double> orientation)
            : D(D), Drot(Drot), position(position), orientation(orientation) {
        type = 0;
        orientvector = vec3<double>(0., 0., 1.);
        orientvector = msmrdtools::rotateVec(orientvector, orientation);
        nextPosition = 1.0 * position;
        nextOrientation = 1.0 * orientation;
        nextOrientvector = 1.0 * orientvector;
    };

    particle::particle(double D, double Drot, std::vector<double> &position, std::vector<double> &orientation)
            : D(D), Drot(Drot), position(position), orientation(orientation) {
        type = 0;
        orientvector = vec3<double>(0., 0., 1.);
        orientvector = msmrdtools::rotateVec(orientvector, orientation);
        nextPosition = vec3<double>(position);
        nextOrientation = quaternion<double>(orientation);
        nextOrientvector = 1.0 * orientvector;
    };


    /* Constructors of particles with Markovian switch: receive input from vec3/quaternion or std::vector and
     * numpy arrays (through pybind) */
    particle::particle(int type0, int state, double D, double Drot, vec3<double> position,
                       quaternion<double> orientation)
            : state(state), D(D), Drot(Drot), position(position), orientation(orientation) {
        type = type0;
    };

    particle::particle(int type0, int state, double D, double Drot, std::vector<double> &position,
                       std::vector<double> &orientation)
            : state(state), D(D), Drot(Drot), position(position), orientation(orientation) {
        type = type0;
    };


    // Implementation of normal particles

    void particle::updatePosition() {
        position = 1.0 * nextPosition;
    };

    void particle::updateOrientation() {
        orientvector = 1 * nextOrientvector;
        orientation = 1 * nextOrientation;
    };

    void particle::setOrientationPyBind(std::vector<double> neworientation) {
        quaternion<double> quat(neworientation);
        orientation = quat;
    }



    //Implementation of particle class with Markovian switch


    // Sets unbound state of particle (the state corresponding to its independent conformation)
    void particle::setState(int newstate) {
        state = newstate;
        boundState = -1;
        boundTo = -1;
    }

    // Sets bound state of particle
    void particle::setBoundState(int newBoundState) {
        state = -1;
        boundState = newBoundState;
        activeMSM = false;
    }



}
