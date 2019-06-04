//
// Created by maojrs on 6/4/19.
//

#include <catch2/catch.hpp>
#include "integrators/integrator.hpp"
#include "integrators/msmrdIntegrator.hpp"
#include "particle.hpp"
#include "quaternion.hpp"
#include "randomgen.hpp"
#include "tools.hpp"
#include "vec3.hpp"

using namespace msmrd;
using msm = msmrd::discreteTimeMarkovStateModel;
using ctmsm = msmrd::continuousTimeMarkovStateModel;
using msmrdMSM = msmrd::msmrdMarkovStateModel;


TEST_CASE("Main MSMRD integrator class", "[msmrdIntegrator]") {
    // Create an msmrdMarkovModel for coupling (tested on testMarkovModels)
    int nBoundStates = 2;
    int nTransitionStates = 2;
    std::map<std::string, float> rateDictionary = { {"1->b1", 5.0}, {"1->b2", 3.0},
                                                    {"2->b1", 4.0}, {"2->b2", 12.0} };
    auto couplingMSM = msmrdMSM(nBoundStates, nTransitionStates, -1, rateDictionary);

    // Create another MSM to model the unbound state dynamics (none in this case)
    int msmid = 0;
    std::vector<std::vector<double>> tmatrix ={ {0} };
    ctmsm unboundMSM = ctmsm(msmid, tmatrix, 0);

    // Create an MSMRD integrator
    double dt = 0.05;
    std::string bodytype = "rigidbody";
    int numParticleTypes = 1;
    double relativeDistanceCutOff = 2.2;
    msmrdIntegrator<ctmsm> testIntegrator = msmrdIntegrator<ctmsm>(dt, 0, bodytype, numParticleTypes,
            relativeDistanceCutOff, unboundMSM, couplingMSM);

    // Test rate dictionary is correctly loaded from couplingMSM into integrator
    REQUIRE(testIntegrator.getRateFromKey("1->b1") == 5.0);
    REQUIRE(testIntegrator.getRateFromKey("2->b1") == 4.0);
    REQUIRE(testIntegrator.getRateFromKey("1->b2") == 3.0);
    REQUIRE(testIntegrator.getRateFromKey("2->b2") == 12.0);

    // Test partition is loaded correctly into integrator
    auto numSections = testIntegrator.positionOrientationPart->getNumSections();
    auto totalSections = numSections[0];
    REQUIRE(totalSections == 203);
    // TEST MORE FUNCTIONS FROM PARTITION AND ALSO IF THE PARTITON IS SET TO A NEW VALUE

    // Tests transition to bound state is correctly computed
    auto transition = testIntegrator.markovModel.computeTransition2BoundState(2);
    auto time = std::get<0>(transition);
    auto endState = std::get<1>(transition);
    REQUIRE(time > 0);
    bool correctEndState = endState == 1 || endState == 2;
    REQUIRE(correctEndState);


}