//
// Created by maojrs on 8/16/18.
//

#include "integrators/overdampedLangevin.hpp"
#include "tools.hpp"

namespace msmrd {
    /**
     * Implementation of over-damped Lanegvin dynamics integrator class
     * @param dt time step
     * @param seed random generator seed (Note seed = -1 corresponds to random device)
     * @param rotation boolean to indicate if rotational degrees of freedom should be integrated
     */
    overdampedLangevin::overdampedLangevin(double dt, long seed, std::string particlesbodytype)
            : integrator(dt, seed, particlesbodytype) {};


    // Integrate one particle from the particle list main routine (visible only inside the class)
    void overdampedLangevin::integrateOne(int partIndex, std::vector<particle> &parts, double timestep) {
        vec3<double> force;
        vec3<double> torque;
        force = 1.0*forceField[partIndex];
        torque = 1.0*torqueField[partIndex];
        translate(parts[partIndex], force, timestep);
        if (rotation) {
            rotate(parts[partIndex], torque, timestep);
        }
    }

    void overdampedLangevin::translate(particle &part, vec3<double> force, double dt0) {
        vec3<double> dr;
        dr = force * dt0 * part.D / KbTemp + std::sqrt(2 * dt0 * part.D) * randg.normal3D(0, 1);
        part.setNextPosition(part.position + dr);
    }

    void overdampedLangevin::rotate(particle &part, vec3<double> torque, double dt0) {
        vec3<double> dphi;
        quaternion<double> dquat;
        dphi = torque * dt0 * part.Drot / KbTemp + std::sqrt(2 * dt0 * part.Drot) * randg.normal3D(0, 1);
        dquat = msmrdtools::axisangle2quaternion(dphi);
        part.setNextOrientation(dquat * part.orientation);
        // Updated orientation vector, useful with rodlike particles
        if (particlesbodytype == "rod" || particlesbodytype == "rodMix") {
            vec3<double> neworientvector = msmrdtools::rotateVec(part.orientvector, dquat);
            part.setNextOrientVector(neworientvector);
        }
    }

}
