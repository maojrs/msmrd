//
// Created by maojrs on 11/25/21.
//

#include "trajectories/trajectoryEnergyTemperature.hpp"

namespace msmrd {

    /**
     * Implementation of trajectory class to store full position only trajectories
     */
    trajectoryEnergyTemperature::trajectoryEnergyTemperature(int bufferSize) : trajectory(1, bufferSize) {
        trajectoryData.reserve(bufferSize);
    }

    // Sample from list of particles and store in trajectoryData
    void trajectoryEnergyTemperature::sample(double time, std::vector<particle> &particleList) {
        std::vector<double> sample(5);
        double kineticEnergy = 0;
        double potentialEnergy = 0;
        double instantTemperature = 0;
        // Calculate kinetic energy
        for (int i = 0; i < particleList.size(); i++) {
            kineticEnergy += 0.5 * particleList[i].mass * particleList[i].velocity.normSquared()/2.0;
        }

        // Calculate instant temperature (assuming Kb=1)
        instantTemperature = kineticEnergy/(3*particleList.size() - 3);

        // Calculate potential energy
        if (externalPotentialActive) {
            for (int i = 0; i < particleList.size(); i++) {
                potentialEnergy += externalPot->evaluate(particleList[i]);;
            }
        }
        if (pairPotentialActive) {
            for (int i = 0; i < particleList.size(); i++) {
                for (int j = i+1; j < particleList.size(); j++) {
                    potentialEnergy += pairPot->evaluate(particleList[i], particleList[j]);
                }
            }
        }
        sample[0] = time;
        sample[1] = kineticEnergy;
        sample[2] = potentialEnergy;
        sample[3] = kineticEnergy + potentialEnergy;
        sample[4] = instantTemperature;
        trajectoryData.push_back(sample);
    }

}
