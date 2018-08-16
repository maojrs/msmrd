//
// Created by maojrs on 8/16/18.
//

#include "integrators/odLangevin.hpp"

/**
 * Class implementation of over-damped Lanegvin dynamics integrator
 * @param dt time step
 * @param seed random generator seed (Note seed = -1 corresponds to random device)
 * @param rotation boolean to indicate if rotational degrees of freedom should be integrated
 */
odLangevin::odLangevin(double dt, long seed, bool rotation) : integrator(dt,seed, rotation) {
    nullPotential nullpot;
    potential0 = &nullpot;
};

void odLangevin::setPotential(potential *pot) {
    potential0 = pot;
}

//void odLangevin::getPotential(std::vector<double> pos) {
//    return potential.evaluate(pos);
//}

// One particle integrate main routine (visible only inside the class)
void odLangevin::integrateOne(particle &part) {
    translate(part,dt);
    if (rotation) {
        rotate(part, dt);
    }
}

// Same as integrateOne, but updates time and publicly visible
void odLangevin::integrate(particle &part) {
    integrateOne(part);
    clock += dt;
}


void odLangevin::translate(particle &part, double dt0){
    vec3<double> dr;
    dr = std::sqrt(2*dt0*part.D)*randg.normal3D(0,1);
    part.setPosition(part.position + dr);
}

void odLangevin::rotate(particle &part, double dt0){
    vec3<double> dphi;
    quaternion<double> dquat;
    dphi = std::sqrt(2*dt0*part.Drot)*randg.normal3D(0,1);
    dquat = angle2quaternion(dphi);
    part.setOrientation(dquat * part.orientation);
}
