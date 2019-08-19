//
// Created by maojrs on 2/6/19.
//

#pragma once
#include "markovModels/markovModel.hpp"

namespace msmrd {
    /**
    * Continuous time Markov state model (ctmsm) class declaration
    */
    class continuousTimeMarkovStateModel : public markovModel {
    private:
        std::vector<double> lambda0;
        std::vector<std::vector<double>> ratescumsum;

        void calculateParameters();

    public:
        continuousTimeMarkovStateModel(int msmid, std::vector<std::vector<double>> tmatrix, long seed);

        void propagate(particleMS &part, int ksteps) override;

        std::tuple<double, int> propagateMSM(int initialState, int ksteps);

        void propagateNoUpdate(particleMS &part, int ksteps);

        // Getters for testing purposes
        std::vector<double> getLambda0() const { return lambda0; };

        std::vector<std::vector<double>> getRatescumsum() const { return ratescumsum; };

    };
}
