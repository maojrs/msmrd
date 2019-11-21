//
// Created by maojrs on 6/4/19.
//

#include <catch2/catch.hpp>
#include "integrators/integrator.hpp"
#include "integrators/msmrdIntegrator.hpp"
#include "boundaries/box.hpp"
#include "discretizations/positionOrientationPartition.hpp"
#include "markovModels/msmrdMarkovModel.hpp"
#include "particle.hpp"
#include "quaternion.hpp"
#include "randomgen.hpp"
#include "tools.hpp"
#include "vec3.hpp"

using namespace msmrd;
using msm = msmrd::discreteTimeMarkovStateModel;
using ctmsm = msmrd::continuousTimeMarkovStateModel;
using msmrdMSM = msmrd::msmrdMarkovModel;

TEST_CASE("Main MSMRD integrator class", "[msmrdIntegrator]") {
    int numBoundStates = 2;
    int numTransitionStates = 2;
    int maxNumBoundStates = 10;
    std::array<double,2> radialBounds{1.25, 2.25};
    double relativeDistanceCutOff = radialBounds[1];
    int numParticleTypes = 1;
    double boxsize = 6;

    // Create MSM for unbound dynamics (no independent conformation switching)
    std::vector<std::vector<double>> tmatrix = {{0}};
    int msmid = 0;
    long seed = -1;
    ctmsm unboundMSM = ctmsm(msmid, tmatrix, seed);
    std::vector<double> Dlist{1.0};
    std::vector<double> Drotlist{1.0};
    unboundMSM.setD(Dlist);
    unboundMSM.setDrot(Drotlist);

    // Create msmrd MSM
    double lagtime = 1;
    std::vector<std::vector<double>> msmrdTmatrix = {{0.0, 0.3, 0.2, 0.5},
                                                     {0.4, 0.3, 0.1, 0.2},
                                                     {0.1, 0.1, 0.6, 0.2},
                                                     {0.4, 0.2, 0.3, 0.1}};
    std::vector<int> activeSet = {1, 2, 11, 12};
    auto msmrdMSM = msmrdMarkovModel(numBoundStates, maxNumBoundStates, msmrdTmatrix, activeSet, lagtime, seed);

    // Create particle list
    int type = 0;
    int state = 0;
    auto position1 = vec3<double> {0.1, 0.1, 0.1};
    auto position2 = vec3<double> {-0.1, -0.1, -0.1};
    auto orientation = quaternion<double> {1.0, 0.0, 0.0, 0.0};
    particle part1 = particle(type, state, Dlist[state], Drotlist[state], position1, orientation);
    particle part2 = particle(type, state, Dlist[state], Drotlist[state], position2, orientation);
    auto plist = std::vector<particle>{part1,part2};

    // Define boundary
    auto boundary = box(boxsize, boxsize, boxsize, "periodic");

    // Define and set up integrators
    double dt = 0.002;
    std::string bodytype = "rigidbody";
    auto myIntegrator = msmrdIntegrator<ctmsm>(dt, seed, bodytype, numParticleTypes, radialBounds, unboundMSM, msmrdMSM);
    myIntegrator.setBoundary(&boundary);

    // Check  msmrdMSM is loaded correctly
    REQUIRE(myIntegrator.markovModel.getTmatrix() == msmrdTmatrix);
    REQUIRE(myIntegrator.markovModel.getNumberBoundStates() == numBoundStates);
    REQUIRE(myIntegrator.markovModel.activeSet == activeSet);

    // Check unboundMSM is loaded correctly
    REQUIRE(myIntegrator.MSMlist[0].getTmatrix() == tmatrix);
    REQUIRE(myIntegrator.MSMlist[0].getNstates() == 1);
    REQUIRE(myIntegrator.MSMlist[0].Dlist == Dlist);
    REQUIRE(myIntegrator.MSMlist[0].Drotlist == Drotlist);

    // Check default discretization is loaded correctly
    REQUIRE(myIntegrator.positionOrientationPart->numTotalSections == 203);

    // Check discretization reassignment works
    int numSphericalSectionsPos = 6;
    int numRadialSectionsQuat = 4;
    int numSphericalSectionsQuat = 6;
    auto discretization = std::make_shared<positionOrientationPartition> (
            positionOrientationPartition(relativeDistanceCutOff, numSphericalSectionsPos,
                                         numRadialSectionsQuat, numSphericalSectionsQuat));
    discretization->setThetasOffset(M_PI/4.0);
    myIntegrator.setDiscretization(discretization);
    REQUIRE(myIntegrator.positionOrientationPart->numTotalSections == 114);
}