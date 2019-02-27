//
// Created by maojrs on 2/6/19.
//

#pragma once
#include "markovModels/markovModel.hpp"

namespace msmrd {
    /**
     * Discrete-time Markov state model class declaration
     */
    class discreteTimeMarkovStateModel : public markovModel {
    public:
        discreteTimeMarkovStateModel(int msmid, std::vector <std::vector<double>> &tempmatrix, double lagtime,
                                     long seed);

        void propagate(particleMS &part, int ksteps) override;

        void propagateNoUpdate(particleMS &part, int ksteps);


    };

}