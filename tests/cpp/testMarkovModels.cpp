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
std::vector<std::vector<double>> tmatrix ={ {-5, 2, 3,}, {3, -6, 3}, {1, 3, -4} };
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
particleMS partMS = particleMS(type, state, D, Drot, position, orientation);
// Create a second particle and a second ctmsm with same seed
ctmsm ctmsmTest2 = ctmsm(msmid, tmatrix, seed);
particleMS partMS2 = particleMS(type, state, D, Drot, position, orientation);
/* Propagate each particle using the same seed but the update and noupdate method,
 * respectively, and compare output at each timestep. */
for (int i=0; i<100; i++) {
ctmsmTest.propagate(partMS, 3);
ctmsmTest2.propagateNoUpdate(partMS2, 3);
partMS2.updateState();
REQUIRE(partMS.getState() == partMS2.getState());
}
}

TEST_CASE("Initialization of ctmsm dictionary in msmrdMarkovModel class", "[msmrdMarkovModel]") {
int nstates = 3;
int numDiscreteOrientations = 6;
std::map<std::string, float> rateDictionary = { {"11->b1", 5.0}, {"11->b2", 3.0},
                                                {"12->b1", 4.0}, {"12->b2", 12.0} };
std::vector<std::vector<double>> tmatrix1 = { {-8, 5, 3}, {0, 0, 0}, {0, 0, 0} };
std::vector<std::vector<double>> tmatrix2 = { {-16, 4, 12}, {0, 0, 0}, {0, 0, 0} };
// Create an msmMarkovModel
auto myMSM = msmrdMarkovModel(nstates, numDiscreteOrientations, -1, rateDictionary);
// Check if the initialization of ctmsms in msmrdMarkovmodel is done correctly.
auto extractedMSM1 = myMSM.getMSM("11");
auto extractedMSM2 = myMSM.getMSM("12");
auto extractedTmatrix1 = extractedMSM1.getTmatrix();
auto extractedTmatrix2= extractedMSM2.getTmatrix();
REQUIRE(extractedTmatrix1 == tmatrix1);
REQUIRE(extractedTmatrix2 == tmatrix2);
}