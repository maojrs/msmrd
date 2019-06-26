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
    double effectiveLagtime = 1.0;
    std::map<std::string, float> rateDictionary = { {"1->b1", 5.0}, {"1->b2", 3.0},
                                                    {"2->b1", 4.0}, {"2->b2", 12.0} };
    auto couplingMSM = msmrdMSM(effectiveLagtime, nBoundStates, nTransitionStates, -1, rateDictionary);

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
    auto relPosition = vec3<double> {0.3, 0.5, 0.7};
    auto relOrientation = quaternion<double> {0.2, 0.4, 0.3, 0.7};
    relOrientation = relOrientation/relOrientation.norm();
    auto refOrientation = quaternion<double> {1.0, 0.0 ,0.0, 0.0};
    auto transitionState = testIntegrator.positionOrientationPart->getSectionNumber(relPosition,
            relOrientation, refOrientation);
    REQUIRE(transitionState == 23);

    // Test if partition works when a new partition value is set
    int numSphericalSectionsPos = 6;
    int numRadialSectionsQuat = 4;
    int numSphericalSectionsQuat = 5;
    positionOrientationPartition newPartition = positionOrientationPartition(relativeDistanceCutOff, numSphericalSectionsPos,
                                                     numRadialSectionsQuat, numSphericalSectionsQuat);
    testIntegrator.setDiscretization(&newPartition);
    numSections = testIntegrator.positionOrientationPart->getNumSections();
    totalSections = numSections[0];
    REQUIRE(totalSections == 96); // total sections =  numSphSecsPos*(numSphSecsQuat*(numRadSecsQuat -1) + 1)


    // Tests transition to bound state is correctly computed by markovModel
    auto transition = testIntegrator.markovModel.computeTransition2BoundState(2);
    auto time = std::get<0>(transition);
    auto endState = std::get<1>(transition);
    REQUIRE(time > 0);
    bool correctEndState = (endState == 1 || endState == 2);
    REQUIRE(correctEndState);

    // Tests transitions to bound states are correctly computed by Integrator
    int type = 0;
    int state = 1;
    double D = 1.0;
    double Drot = 1.0;
    auto position1 = vec3<double> {0.1, 0.1, 0.1};
    auto position2 = vec3<double> {-0.1, -0.1, -0.1};
    auto orientation = quaternion<double> {1.0, 0.0, 0.0, 0.0};
    particleMS partMS1 = particleMS(type, state, D, Drot, position1, orientation);
    particleMS partMS2 = particleMS(type, state, D, Drot, position2, orientation);
    auto plist = std::vector<particleMS>{partMS1,partMS2};
    auto numEvents = testIntegrator.eventMgr.getNumEvents();
    REQUIRE(numEvents == 0);
    testIntegrator.computeTransitions2BoundStates(plist);
    numEvents = testIntegrator.eventMgr.getNumEvents();
    auto event = testIntegrator.eventMgr.getEvent(0,1);
    auto transitionTime = event.waitTime;
    auto nextState = event.endState;
    auto iIndex = event.part1Index;
    auto jIndex = event.part2Index;
    auto inORout = event.inORout;
    REQUIRE(numEvents == 1);
    REQUIRE(transitionTime > 0);
    REQUIRE( (nextState == 1 || nextState == 2) );
    REQUIRE(iIndex == 0);
    REQUIRE(jIndex == 1);
    REQUIRE( (inORout == "in" || inORout == "out") );

}