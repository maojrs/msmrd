//
// Created by maojrs on 3/28/19.
//
#include <catch2/catch.hpp>
#include "markovModels/discreteTimeMarkovModel.hpp"
#include "markovModels/continuousTimeMarkovModel.hpp"
#include "markovModels/msmrdMarkovModel.hpp"


using namespace msmrd;
using msm = msmrd::discreteTimeMarkovStateModel;
using ctmsm = msmrd::continuousTimeMarkovStateModel;

TEST_CASE("Fundamental CTMSM parameters and propagation test", "[ctmsm]") {
    int msmid = 0;
    std::vector<std::vector<double>> tmatrix ={ {-5, 2, 3}, {3, -6, 3}, {1, 3, -4} };
    // sum of outgoing rates for each state (row sum without diagonal negative value)
    std::vector<double> lambda0 = {5, 6, 4};
    // cumulative sum of outgoing rates for each state
    std::vector<std::vector<double>> ratescumsum = {{2, 5}, {3, 6}, {1, 4}};
    long seed = 0;
    ctmsm ctmsmTest = ctmsm(msmid, tmatrix, seed);
    REQUIRE(ctmsmTest.nstates == 3);
    REQUIRE(ctmsmTest.getLambda0() == lambda0);
    REQUIRE(ctmsmTest.getRatescumsum() == ratescumsum);
    // Propagation test requires creating a particleMS
    int type = 0;
    int state = 2;
    double D = 1.0;
    double Drot = 0.5;
    auto position = vec3<double> {0.0, 0.0, 0.0};
    auto orientation = quaternion<double> {1.0, 0.0, 0.0, 0.0};
    particle partMS = particle(type, state, D, Drot, position, orientation);
    // Create a second particle and a second ctmsm with same seed
    ctmsm ctmsmTest2 = ctmsm(msmid, tmatrix, seed);
    particle partMS2 = particle(type, state, D, Drot, position, orientation);
    /* Propagate each particle using the same seed but the update and noupdate method,
     * respectively, and compare output at each timestep. */
    for (int i=0; i<100; i++) {
        ctmsmTest.propagate(partMS, 3);
        ctmsmTest2.propagateNoUpdate(partMS2, 3);
        partMS2.updateState();
        REQUIRE(partMS.getState() == partMS2.getState());
    }
}

TEST_CASE("Initialization of msmrdMarkovModel class", "[msmrdMarkovModel]") {
    double lagtime = 1;
    long seed = -1;
    int numBoundStates = 2; // states 1 and 2
    int maxNumBoundStates = 10;
    int numTransitionStates = 2; // states 10 or 11
    std::vector<std::vector<double>> tmatrix = {{0.0, 0.3, 0.2, 0.5},
                                                {0.4, 0.3, 0.1, 0.2},
                                                {0.1, 0.1, 0.6, 0.2},
                                                {0.4, 0.2, 0.3, 0.1}};
    std::vector<int> activeSet = {1, 2, 11, 12}; // MSMindexing is 0, 1, 2, 3, respectively.

    auto msmrdMSM = msmrdMarkovModel(numBoundStates, maxNumBoundStates, tmatrix, activeSet, lagtime, seed);

    // Check transition works normally
    auto transition = msmrdMSM.calculateTransition(1); // uses activeSetIndex
    auto time = std::get<0>(transition);
    auto endState = std::get<1>(transition); // gets state in activeSetIndex
    REQUIRE(time > 0);
    bool correctEndState = endState == 2 || endState == 11 || endState == 12; // uses activeSetIndex
    REQUIRE(correctEndState);
    // Check functions to recover MSM and activeSet indexes are working
    for (int i = 0; i < activeSet.size(); i++) {
        REQUIRE(msmrdMSM.getActiveSetIndex(i) == activeSet[i]);
        REQUIRE(msmrdMSM.getMSMindex(activeSet[i]) == i);
    }
}