//
// Created by maojrs on 4/6/22.
//
#include "potentials/combinedPotentials.hpp"

namespace msmrd {
    /*
     * Implementations for combined potentials
     */

    void combinedPairPotential::addPotential(pairPotential *pairPotentialPtr) {
        potentials.push_back(std::move(std::shared_ptr<pairPotential>(pairPotentialPtr)));
    }

    double combinedPairPotential::evaluate(particle &part1, particle &part2) {
        double potentialVal = 0.0;
        for(const auto& potential: potentials) {
            potentialVal += potential->evaluate(part1, part2);
        }
        return potentialVal;
    };

    std::array<vec3<double>, 4> combinedPairPotential::forceTorque(particle &part1, particle &part2) {
        vec3<double> force = vec3<double>(0, 0, 0);
        vec3<double> torque = vec3<double>(0, 0, 0);
        for(const auto& potential: potentials) {
            auto forctorq = potential->forceTorque(part1, part2);
            force += forctorq[0];
            torque += forctorq[1];
        }
        return {force, torque, -1.0*force, -1.0*torque};
    };


    void combinedExternalPotential::addPotential(externalPotential *externalPotentialPtr) {
        potentials.push_back(std::move(std::shared_ptr<externalPotential>(externalPotentialPtr)));
        externalPotentialPtr;
    }

    double combinedExternalPotential::evaluate(particle &part1) {
        double potentialVal = 0.0;
        for(const auto& potential: potentials) {
            potentialVal += potential->evaluate(part1);
        }
        return potentialVal;
    };

    std::array<vec3<double>, 2> combinedExternalPotential::forceTorque(particle &part1) {
        vec3<double> force = vec3<double>(0, 0, 0);
        vec3<double> torque = vec3<double>(0, 0, 0);
        for(const auto& potential: potentials) {
            auto forctorq = potential->forceTorque(part1);
            force += forctorq[0];
            torque += forctorq[1];
        }
        return {force, torque};
    };

}