//
// Created by maojrs on 4/19/21.
//

#include "potentials/lennardJones.hpp"

namespace msmrd {
    /**
     * Class implementation for Lennard Jonnes pair potential
     */
    lennardJones::lennardJones(double epsilon, double sigma) : epsilon(epsilon), sigma(sigma) {
        cutOff = 3 * sigma;
    }

    lennardJones::lennardJones(double epsilon, double sigma, double cutOff) : epsilon(epsilon), sigma(sigma),
        cutOff(cutOff) {}

    lennardJones::lennardJones(double epsilon, double sigma,  std::vector<int> exclParticleTypesPairs)
            : epsilon(epsilon), sigma(sigma)    {
        cutOff = 3 * sigma;
        excludeParticleTypesPairs = exclParticleTypesPairs;
    }

    /* Set maximum value of force. If potnetial yields larger force, uses this value. Also calculates
     * the corresponding value at which the potential is cut off. */
    void lennardJones::setForceCapValue(double forceCapVal) {
        forceCap = false;
        forceCapValue = forceCapVal;
        // Calculate potential cap
        double rr =  cutOff;
        double drr = rr/1000000.0;
        double forceNorm = 0.0;
        auto position = vec3<double>{0.0, 0.0, 0.0};
        auto orientation = quaternion<double>{1.0, 0.0, 0.0, 0.0};
        particle part1 = particle(0, 0, 0, 0, position, orientation);
        particle part2 = particle(0, 0, 0, 0, position, orientation);
        while (forceNorm <= forceCapValue) {
            part2.setPosition(vec3<double>{rr, 0, 0});
            auto force = forceTorque(part1, part2);
            forceNorm = force[0].norm();
            rr = rr - drr;
            if (rr < 0) {
                break;
            }
        }
        rcritical = 1.0 * rr;
        potentialCutOff = 1.0 * evaluate(part1, part2);
        forceCap = true;
    }

    /* Set a value at which the repulsive potential is cutoff and substituted by a line with constant
     * slope (correponding to constant force). Also sets the corresponding force cap value. */
    void lennardJones::setPotentialCutOff(double potCutOff) {
        forceCap = false;
        potentialCutOff = potCutOff;
        // Calculate force cap (solve quadratic)
        double termInside = (2 * epsilon / potentialCutOff) * (-1 + std::sqrt(1 + potentialCutOff/epsilon));
        rcritical = sigma * std::pow(termInside, 1.0/6.0);
        auto position1 = vec3<double>{0.0, 0.0, 0.0};
        auto position2 = vec3<double>{rcritical, 0.0, 0.0};
        auto orientation = quaternion<double>{1, 0.0, 0.0, 0.0};
        particle part1 = particle(0, 0, 0, 0, position1, orientation);
        particle part2 = particle(0, 0, 0, 0, position2, orientation);
        auto force = forceTorque(part1, part2);
        forceCapValue = force[0].norm();
        forceCap = true;
    }

    double lennardJones::evaluate(particle &part1, particle &part2) {
        std::array<vec3<double>, 2> relPos = relativePositionComplete(part1.position, part2.position);
        //vec3<double> pos1virtual = relPos[0]; // virtual pos1 if periodic boundary; otherwise pos1.
        vec3<double> rvec = relPos[1]; //part2.position - part1.position if non-periodic boundary;
        double r = rvec.norm();
        double resultingPotential = 0;
        auto activeParticles = true;
        if (r <= cutOff) {
            if (std::find(excludeParticleTypesPairs.begin(),
                          excludeParticleTypesPairs.end(), part1.type) == excludeParticleTypesPairs.end() and
                std::find(excludeParticleTypesPairs.begin(),
                          excludeParticleTypesPairs.end(), part2.type) == excludeParticleTypesPairs.end()) {
                activeParticles = false;
            }
            if (activeParticles) {
                auto term0 = std::pow((sigma / r), 6);
                auto term1 = std::pow(term0, 2);
                resultingPotential = 4 * epsilon * (term1 - term0) + baseEnergy;
                if (forceCap and resultingPotential >= potentialCutOff) {
                    resultingPotential = -forceCapValue * (r - rcritical) + potentialCutOff;
                }
            }
        }
        return resultingPotential;
    }

    std::array<vec3<double>, 4> lennardJones::forceTorque(particle &part1, particle &part2) {
        std::array<vec3<double>, 2> relPos = relativePositionComplete(part1.position, part2.position);
        //vec3<double> pos1virtual = relPos[0]; // virtual pos1 if periodic boundary; otherwise pos1.
        vec3<double> rvec = relPos[1]; //part2.position - part1.position if non-periodic boundary;
        double r = rvec.norm();
        vec3<double> force = vec3<double>();
        vec3<double> torque = vec3<double>();
        auto activeParticles = true;
        if (r <= cutOff) {
            if (std::find(excludeParticleTypesPairs.begin(),
                          excludeParticleTypesPairs.end(), part1.type) == excludeParticleTypesPairs.end() and
                std::find(excludeParticleTypesPairs.begin(),
                          excludeParticleTypesPairs.end(), part2.type) == excludeParticleTypesPairs.end()) {
                activeParticles = false;
            }
            if (activeParticles) {
                double term0 = std::pow((sigma / r), 6);
                double dVdr = (24 * epsilon / r) * (term0 - 2 * std::pow(term0, 2));
                double dVdx = dVdr * rvec[0] / r;
                double dVdy = dVdr * rvec[1] / r;
                double dVdz = dVdr * rvec[2] / r;

                force = vec3<double>(dVdx, dVdy, dVdz);
                torque = 0 * force;
                // Apply force cap if active
                if (forceCap and force.norm() >= forceCapValue) {
                    force = forceCapValue * force / force.norm();
                }
                return {force, torque, -1 * force, torque};
            }
        }
        return {force, torque, -1 * force, torque};
    }


    /**
     * Class implementation for WCA pair potential
     */
    WCA::WCA(double epsilon, double sigma) : lennardJones(epsilon,sigma) {
        cutOff = std::pow(2,1.0/6.0) * sigma;
        baseEnergy = 1 * epsilon;
    }

    WCA::WCA(double epsilon, double sigma, std::vector<int> exclParticleTypesPairs)
            :lennardJones(epsilon,sigma,exclParticleTypesPairs) {
        cutOff = std::pow(2,1.0/6.0) * sigma;
        baseEnergy = 1 * epsilon;
    }


}