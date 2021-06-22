//
// Created by maojrs on 6/14/21.
//

#include "trajectories/trajectoryMoriZwanzig.hpp"

namespace msmrd {

    void trajectoryMoriZwanzig::sample(double time, std::vector<particle> &particleList) {
        std::vector<double> sample(8);
        for (int i = 0; i < particleList.size(); i++) {
            // If particle type corresponds to one of the distinguished particles, sample its value
            if (std::find(distinguishedTypes.begin(), distinguishedTypes.end(),
                          particleList[i].type) != distinguishedTypes.end()) {
                sample[0] = time;
                for (int k = 0; k < 3; k++) {
                    sample[k + 1] = particleList[i].position[k];
                }
                sample[4] = 1.0 * particleList[i].type;
                for (int k = 0; k < 3; k++) {
                    sample[k + 5] = particleList[i].raux[k];
                }
                trajectoryData.push_back(sample);
            };
        }
    }

}