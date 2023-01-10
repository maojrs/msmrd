//
// Created by maojrs on 2/8/22.
//

#include "vec3.hpp"
#include "potentials/pairBistable.hpp"

namespace msmrd {
    /* Class implementation for pair bistable potential.
     * The potential is of the form: scaleFactor ( 1 - ((x-(x0+rad))/rad)^2 )^2 */

    pairBistable::pairBistable(double x0, double rad, double scalefactor) : x0(x0), rad(rad),
    scalefactor(scalefactor) {}

    pairBistable::pairBistable(double x0, double rad, std::vector<int> partTypes,
            double scalefactor) : pairBistable(x0, rad, scalefactor) {
        particleTypes = partTypes;
    }

    // Evaluate potential value for two given particles' positions
    double pairBistable::evaluate(particle &part1, particle &part2) {
        bool activeParticles = false;
        double potential = 0.0;
        if (std::find(particleTypes.begin(),particleTypes.end(), part1.type) != particleTypes.end() and
            std::find(particleTypes.begin(),particleTypes.end(), part2.type) != particleTypes.end()) {
            activeParticles = true;
        }
        if (particleTypes.empty() or activeParticles) {
            vec3<double> rvec = relativePosition(part2.position, part1.position); //part1.position - part2.position;
            double arg = (rvec.norm() - (x0 + rad))/rad;
            potential =  scalefactor * std::pow(1 - std::pow(arg,2) ,2);
        }
        return potential;
    }

    // Returns -gradient of potential (force)  and zero torque at position x
    std::array<vec3<double>, 4> pairBistable::forceTorque(particle &part1, particle &part2) {
        vec3<double> force = vec3<double>(0, 0, 0);
        vec3<double> torque = vec3<double>(0, 0, 0);
        bool activeParticles = false;
        if (std::find(particleTypes.begin(),particleTypes.end(), part1.type) != particleTypes.end() and
            std::find(particleTypes.begin(),particleTypes.end(), part2.type) != particleTypes.end()) {
            activeParticles = true;
        }
        if (particleTypes.empty() or activeParticles) {
            vec3<double> rvec = relativePosition(part2.position, part1.position); //part1.position - part2.position;
            double r = rvec.norm();
            double arg = (r - (x0 + rad))/rad;
            double dVdr = 4.0 * (1 - std::pow(arg,2)) * arg /rad;
            double dVdx = dVdr * rvec[0] / r;
            double dVdy = dVdr * rvec[1] / r;
            double dVdz = dVdr * rvec[2] / r;
            force = vec3<double>(dVdx, dVdy, dVdz);
            force = scalefactor * force;
        }
        return {force, torque, -1.0*force, -1.0*torque};
    }


    /* Class implementation for pair bistable potential with a bias.
     * The potential is of the form: scaleFactor ( 1 - ((x-(x0+rad))/rad)^2 )^2 + a*log(x)*/

    pairBistableBias::pairBistableBias(double x0, double rad, double scalefactor) : x0(x0), rad(rad),
                                                                            scalefactor(scalefactor) {}

    pairBistableBias::pairBistableBias(double x0, double rad, std::vector<int> partTypes,
                               double scalefactor) : pairBistableBias(x0, rad, scalefactor) {
        particleTypes = partTypes;
    }

    // Evaluate potential value for two given particles' positions
    double pairBistableBias::evaluate(particle &part1, particle &part2) {
        bool activeParticles = false;
        double potential = 0.0;
        if (std::find(particleTypes.begin(),particleTypes.end(), part1.type) != particleTypes.end() and
            std::find(particleTypes.begin(),particleTypes.end(), part2.type) != particleTypes.end()) {
            activeParticles = true;
        }
        if (particleTypes.empty() or activeParticles) {
            vec3<double> rvec = relativePosition(part2.position, part1.position); //part1.position - part2.position;
            double r = rvec.norm();
            double arg = (r - (x0 + rad))/rad;
            potential =  scalefactor * (std::pow(1 - std::pow(arg,2) ,2) + a * std::log(r) );
        }
        return potential;
    }

    // Returns -gradient of potential (force)  and zero torque at position x
    std::array<vec3<double>, 4> pairBistableBias::forceTorque(particle &part1, particle &part2) {
        vec3<double> force = vec3<double>(0, 0, 0);
        vec3<double> torque = vec3<double>(0, 0, 0);
        bool activeParticles = false;
        if (std::find(particleTypes.begin(),particleTypes.end(), part1.type) != particleTypes.end() and
            std::find(particleTypes.begin(),particleTypes.end(), part2.type) != particleTypes.end()) {
            activeParticles = true;
        }
        if (particleTypes.empty() or activeParticles) {
            vec3<double> rvec = relativePosition(part2.position, part1.position); //part1.position - part2.position;
            double r = rvec.norm();
            double arg = (r - (x0 + rad))/rad;
            double dVdr = 4.0 * (1 - std::pow(arg,2)) * arg /rad + a/r;
            double dVdx = dVdr * rvec[0] / r;
            double dVdy = dVdr * rvec[1] / r;
            double dVdz = dVdr * rvec[2] / r;
            force = vec3<double>(dVdx, dVdy, dVdz);
            force = scalefactor * force;
        }
        return {force, torque, -1.0*force, -1.0*torque};
    }

}