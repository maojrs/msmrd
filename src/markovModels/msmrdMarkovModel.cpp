//
// Created by maojrs on 2/6/19.
//

#include "markovModels/msmrdMarkovModel.hpp"

namespace msmrd{
    using ctmsm = msmrd::continuousTimeMarkovStateModel;
    /**
    * Implementation of Markov model class specializes for MSM/RD.
    * @param rateDict
    */
    msmrdMarkovModel::msmrdMarkovModel(unsigned int nstates, unsigned int numDiscreteOrientations,
                                       long seed, std::map<std::string, float> &rateDictionary)
            : nstates(nstates), numDiscreteOrientations(numDiscreteOrientations),
              seed(seed), rateDictionary(rateDictionary) {
        generateMarkovModels();
    };

    /* Generates several continuous-time Markov models each representing a discretization of
     * the relative orientation */
    void msmrdMarkovModel::generateMarkovModels() {
        int numBoundStates = nstates - 1;
        //int numMSMs = numBoundStates + (numDiscreteOrientations)*(numDiscreteOrientations + 1)/2;

        // Create all MSMs corresponding to different orientations
        float rate;
        std::vector<std::vector<double>> transitionMatrix;
        int msmid = 0;
        for (int i = 1; i <= numDiscreteOrientations; i++){
            for (int j = i; j <= numDiscreteOrientations; j++){
                // Reste values of transitio matrix to zero.
                transitionMatrix = std::vector<std::vector<double>> (nstates, std::vector<double> (nstates, 0));
                for (int k = 1; k <= numBoundStates; k++) {
                    // Extract on rates to bound states
                    auto key = std::to_string(i) + std::to_string(j) + "->b" + std::to_string(k);
                    auto search = rateDictionary.find(key);
                    if (search != rateDictionary.end()) {
                        rate = search->second;
                    } else {
                        rate = 0;
                    }
                    transitionMatrix[0][k] = rate;
                    transitionMatrix[0][0] -= rate;
                }
                // Create ctmsm and insert it into dictionary with key 'ij'.
                ctmsm ctmsm_ij = ctmsm(msmid, transitionMatrix, seed*(msmid+1));
                msmid += 1;
                auto newkey = std::to_string(i) + std::to_string(j);
                ctmsmDictionary.insert ( std::pair<std::string, ctmsm>(newkey, ctmsm_ij) );
            }
        }


//        // Extract off rates to unbound transition states
//        auto key2 = "b" + std::to_string(k) + "->" + std::to_string(i) + std::to_string(j);
//        auto search2 = rateDictionary.find(key2);
//        if (search2 != rateDictionary.end()) {
//            rate = search2->second;
//        } else {
//            rate = 0;
//        }
//        transitionMatrix[k][0] = rate;
//        transitionMatrix[k][0] -= rate;

    }
}
