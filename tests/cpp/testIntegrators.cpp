//
// Created by maojrs on 6/4/19.
//

#include <catch2/catch.hpp>
#include "integrators/msmrdIntegrator.hpp"
#include "integrators/msmrdMultiParticleIntegrator.hpp"
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

TEST_CASE("Initialization of MSMRD integrator class", "[msmrdIntegrator]") {
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
    auto myIntegrator = msmrdIntegrator<ctmsm>(dt, seed, bodytype, numParticleTypes,
                                               radialBounds, unboundMSM, msmrdMSM);
    myIntegrator.setBoundary(&boundary);

    // Check  msmrdMSM is loaded correctly
    REQUIRE(myIntegrator.msmrdMSM.getTmatrix() == msmrdTmatrix);
    REQUIRE(myIntegrator.msmrdMSM.getNumberBoundStates() == numBoundStates);
    REQUIRE(myIntegrator.msmrdMSM.activeSet == activeSet);

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
    double offset = M_PI/4.0;
    discretization->setThetasOffset(offset);
    myIntegrator.setDiscretization(discretization);
    REQUIRE(myIntegrator.positionOrientationPart->numTotalSections == 114);
    auto partition = myIntegrator.positionOrientationPart->sphericalPartition->getPartition();
    auto thetas = std::get<2>(partition);
    auto thetasRef = std::vector<double> {0.0 + offset, 1.5707963267948966 + offset,
                                          3.141592653589793 + offset, 4.71238898038469 + offset};
    REQUIRE(thetas[0] == thetasRef);
}

TEST_CASE("Initialization and functions of MSMRD multi-particle integrator class", "[msmrdMultiParticleIntegrator]") {
    int numBoundStates = 4;
    int numTransitionStates = 2;
    int maxNumBoundStates = 10;
    std::array<double,2> radialBounds{1.25, 2.25};
    double relativeDistanceCutOff = radialBounds[1];
    int numParticleTypes = 1;
    double boxsize = 2.5;

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
    std::vector<std::vector<double>> msmrdTmatrix = {{0.0, 0.3, 0.2, 0.5, 0.0, 0.0},
                                                     {0.4, 0.3, 0.1, 0.2, 0.0, 0.0},
                                                     {0.1, 0.1, 0.6, 0.2, 0.0, 0.0},
                                                     {0.4, 0.2, 0.3, 0.1, 0.0, 0.0},
                                                     {0.2, 0.1, 0.0, 0.0, 0.4, 0.3},
                                                     {0.0, 0.0, 0.1, 0.3, 0.2, 0.4}};
    std::vector<int> activeSet = {1, 2, 11, 12};
    auto msmrdMSM = msmrdMarkovModel(numBoundStates, maxNumBoundStates, msmrdTmatrix, activeSet, lagtime, seed);

    // Create particle list
    int type = 0;
    int state = 0;
    auto position0 = vec3<double> {0.1, 0.1, 0.1};
    auto position1 = vec3<double> {-0.1, -0.1, -0.1};
    auto position2 = vec3<double> {-0.1, 0.1, -0.1};
    auto position3 = vec3<double> {0.1, 0.1, -0.1};
    auto position4 = vec3<double> {-0.1, 0.1, 0.1};
    auto orientation = quaternion<double> {1.0, 0.0, 0.0, 0.0};
    particle part0 = particle(type, state, Dlist[state], Drotlist[state], position0, orientation);
    particle part1 = particle(type, state, Dlist[state], Drotlist[state], position1, orientation);
    particle part2 = particle(type, state, Dlist[state], Drotlist[state], position2, orientation);
    particle part3 = particle(type, state, Dlist[state], Drotlist[state], position3, orientation);
    particle part4 = particle(type, state, Dlist[state], Drotlist[state], position4, orientation);
    auto plist = std::vector<particle>{part0,part1,part2, part3, part4};

    // Define boundary
    auto boundary = box(boxsize, boxsize, boxsize, "periodic");

    // Define and set up integrators
    double dt = 0.0001;
    std::string bodytype = "rigidbody";
    auto myIntegrator = msmrdMultiParticleIntegrator<ctmsm>(dt, seed, bodytype, numParticleTypes,
            radialBounds, unboundMSM, msmrdMSM);
    myIntegrator.setBoundary(&boundary);

    // Check  msmrdMSM is loaded correctly
    REQUIRE(myIntegrator.msmrdMSM.getTmatrix() == msmrdTmatrix);
    REQUIRE(myIntegrator.msmrdMSM.getNumberBoundStates() == numBoundStates);
    REQUIRE(myIntegrator.msmrdMSM.activeSet == activeSet);

    // Check unboundMSM is loaded correctly
    REQUIRE(myIntegrator.MSMlist[0].getTmatrix() == tmatrix);
    REQUIRE(myIntegrator.MSMlist[0].getNstates() == 1);
    REQUIRE(myIntegrator.MSMlist[0].Dlist == Dlist);
    REQUIRE(myIntegrator.MSMlist[0].Drotlist == Drotlist);

    // Check default discretization is loaded correctly
    REQUIRE(myIntegrator.positionOrientationPart->numTotalSections == 203);

    // Check adding a particle compound (note position of binding particle (jIndex) will change)
    REQUIRE(myIntegrator.particleCompounds.size() == 0);
    // Bind particle 0 and 2 with bound state 0
    int iIndex = 0;
    int jIndex = 2;
    int boundState = 0;
    myIntegrator.addCompound(plist, iIndex, jIndex, boundState);
    REQUIRE(myIntegrator.particleCompounds.size() == 1);
    REQUIRE(myIntegrator.particleCompounds[0].compoundSize == 2);
    // Check particles are set into the correct particleCompound index
    REQUIRE(plist[iIndex]. compoundIndex == 0);
    REQUIRE(plist[jIndex]. compoundIndex == 0);
    // Check the newly bound particle is set into the correct relative position and orientation
    auto relPosition = myIntegrator.discreteTrajClass->getRelativePosition(boundState);
    auto relOrientation = myIntegrator.discreteTrajClass->getRelativeOrientation(boundState);
    REQUIRE(plist[iIndex].position + relPosition == plist[jIndex].position);
    REQUIRE(relOrientation * plist[iIndex].orientation == plist[jIndex].orientation);

    // Check binding another particle into exisiting compound: bind particle 1 and 2 with bound state 3
    iIndex = 1;
    jIndex = 2;
    boundState = 3;
    myIntegrator.addCompound(plist, iIndex, jIndex, boundState);
    REQUIRE(myIntegrator.particleCompounds.size() == 1);
    REQUIRE(myIntegrator.particleCompounds[0].compoundSize == 3);
    // Check particles are set into the correct particleCompound index
    REQUIRE(plist[0].compoundIndex == 0);
    REQUIRE(plist[1].compoundIndex == 0);
    REQUIRE(plist[2].compoundIndex == 0);
    relPosition = myIntegrator.discreteTrajClass->getRelativePosition(boundState);
    relOrientation = myIntegrator.discreteTrajClass->getRelativeOrientation(boundState);
    REQUIRE(plist[jIndex].position + relPosition == plist[iIndex].position);
    REQUIRE(relOrientation * plist[jIndex].orientation == plist[iIndex].orientation);

    // Make another independent complex: bind particle 3 and 4 with bound state 2
    iIndex = 3;
    jIndex = 4;
    boundState = 2;
    myIntegrator.addCompound(plist, iIndex, jIndex, boundState);
    REQUIRE(myIntegrator.particleCompounds.size() == 2);
    REQUIRE(myIntegrator.particleCompounds[1].compoundSize == 2);
    // Check particles are set into the correct particleCompound index
    REQUIRE(plist[iIndex].compoundIndex == 1);
    REQUIRE(plist[jIndex].compoundIndex == 1);
    // Check the newly bound particle is set into the correct relative position and orientation
    relPosition = myIntegrator.discreteTrajClass->getRelativePosition(boundState);
    relOrientation = myIntegrator.discreteTrajClass->getRelativeOrientation(boundState);
    REQUIRE(plist[iIndex].position + relPosition == plist[jIndex].position);
    REQUIRE(relOrientation * plist[iIndex].orientation == plist[jIndex].orientation);

    // Diffuse compounds
    int timesteps = 1;
    auto prevPosition = 1.0 * myIntegrator.particleCompounds[0].position;
    auto prevOrientation = 1.0 * myIntegrator.particleCompounds[0].orientation;
    for (int i = 0; i < timesteps; i++) {
        myIntegrator.integrateDiffusionCompounds(plist,0.001);
    }
    //REQUIRE(prevPosition != myIntegrator.particleCompounds[1].position);
    //REQUIRE(prevOrientation != myIntegrator.particleCompounds[1].orientation);
    // Check relative orientation and orientation is maintained after integration
    auto newRelPosition = msmrdtools::rotateVec(relPosition, myIntegrator.particleCompounds[1].orientation);
    REQUIRE((plist[iIndex].position + newRelPosition - plist[jIndex].position).norm() <= 0.000001);
    REQUIRE((relOrientation * plist[iIndex].orientation - plist[jIndex].orientation).norm() <= 0.000001);
}
