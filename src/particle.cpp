//
// Created by maojrs on 10/17/19.
//

#include "particle.hpp"

namespace msmrd {

    /**
     * Implementation of particle class
     */

    // Constructors: receive input from vec3/quaternion or std::vector and numpy arrays (through pybind)
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


    /**
     * Implementation of particleMS class
     */


    particleMS::particleMS(int type0, int state, double D, double Drot, vec3<double> position,
               quaternion<double> orientation)
            : state(state), particle(D, Drot, position, orientation) {
        type = type0;
    };

    particleMS::particleMS(int type0, int state, double D, double Drot, std::vector<double> &position,
               std::vector<double> &orientation)
            : state(state), particle(D, Drot, position, orientation) {
        type = type0;
    };

    // Sets unbound state of particle (the state corresponding to its independent conformation)
    void particleMS::setState(int newstate) {
        state = newstate;
        boundState = -1;
        boundTo = -1;
    }

    // Sets bound state of particle
    void particleMS::setBoundState(int newBoundState) {
        state = -1;
        boundState = newBoundState;
        activeMSM = false;
    }


    /**
     * Implementation of particleComplex class
     */

    /* Consturctors */
    particleCompound::particleCompound(vec3<double> position) : position(position){};

    particleCompound::particleCompound(std::vector<double> &position) : position(position) {};

    particleCompound::particleCompound(std::map<std::tuple<int,int>, int> boundPairsDictionary):
            boundPairsDictionary(boundPairsDictionary) {};

    particleCompound::particleCompound(vec3<double> position, std::map<std::tuple<int,int>, int> boundPairsDictionary):
    position(position), boundPairsDictionary(boundPairsDictionary) {};

    particleCompound::particleCompound(std::vector<double> &position,
                                       std::map<std::tuple<int,int>, int> boundPairsDictionary) :
    position(position), boundPairsDictionary(boundPairsDictionary) {};


    /* Joins another particle complex into this particle complex. The local dictionary has preference if
     * equal keys. However, that should never happen. */
    void particleCompound::joinParticleCompound(particleCompound pComplex) {
        boundPairsDictionary.insert(pComplex.boundPairsDictionary.begin(),
                                    pComplex.boundPairsDictionary.end());
    };



}
