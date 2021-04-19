//
// Created by maojrs on 4/19/21.
//

#include "integrators/langevin.hpp"

namespace msmrd {
    /**
     * Implementation of Lanegvin dynamics integrator class
     * @param dt time step
     * @param seed random generator seed (Note seed = -1 corresponds to random device)
     * @param rotation boolean to indicate if rotational degrees of freedom should be integrated
     */
    langevin::langevin(double dt, long seed, std::string particlesbodytype)
            : integrator(dt, seed, particlesbodytype) {
        rotation = false;
        velocityIntegration = true;
        if (particlesbodytype != "point") {
            throw std::invalid_argument("Langevin integrator only implemented for point particles "
                                        "(no rotation allowed)");
        }
    }

    // Integrate one particle from the particle list main routine (visible only inside the class)
    void langevin::integrateOne(int partIndex, std::vector<particle> &parts, double timestep) {
        vec3<double> force;
        vec3<double> torque;
        force = 1.0*forceField[partIndex];
        torque = 1.0*torqueField[partIndex];
        translate(parts[partIndex], force, timestep);
    }

    // Translate one particle using Euler Maruyama on the Langevin equation
    void langevin::translate(particle &part, vec3<double> force, double dt0) {
        vec3<double> dr;
        vec3<double> dv;
        double eta = KbTemp / part.D; // friction coefficient
        dr = dv * dt0;
        dv = - eta * part.velocity * dt0 / part.mass
                + force * dt0 / part.mass
                + (std::sqrt(2 * KbTemp * eta * dt0 )/part.mass) * randg.normal3D(0, 1);
        part.setNextPosition(part.position + dr);
        part.setNextVelocity(part.velocity + dv);
    }



}