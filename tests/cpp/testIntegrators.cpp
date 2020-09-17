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
    auto orientationA = quaternion<double> {2.0, 1.0, 4.0, 0.5};
    orientationA = orientationA/orientationA.norm();
    auto orientationB = quaternion<double> {0.5, 0.7, 1.0, 2.0};
    orientationB = orientationB/orientationB.norm();
    particle part0 = particle(type, state, Dlist[state], Drotlist[state], position0, orientationA);
    particle part1 = particle(type, state, Dlist[state], Drotlist[state], position1, orientationA);
    particle part2 = particle(type, state, Dlist[state], Drotlist[state], position2, orientationB);
    particle part3 = particle(type, state, Dlist[state], Drotlist[state], position3, orientationB);
    particle part4 = particle(type, state, Dlist[state], Drotlist[state], position4, orientationA);
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
    REQUIRE(myIntegrator.getCompoundSize(0) == 2);
    // Check particles are set into the correct particleCompound index
    REQUIRE(plist[iIndex]. compoundIndex == 0);
    REQUIRE(plist[jIndex]. compoundIndex == 0);
    /* Check the newly bound particle is set into the correct relative position and orientation.
     * Note this relative positions/orientations are w/respect to other particle and not the center
     * of the particle compound.*/
    auto relPosition = myIntegrator.discreteTrajClass->getRelativePosition(boundState);
    auto relOrientation = myIntegrator.discreteTrajClass->getRelativeOrientation(boundState);
    auto newRelPosition = msmrdtools::rotateVec(relPosition, plist[iIndex].orientation);
    REQUIRE(plist[iIndex].position + newRelPosition == plist[jIndex].position);
    REQUIRE(plist[iIndex].orientation * relOrientation == plist[jIndex].orientation);

    // Check relative position of reference particle is correct
    auto relPosRef = plist[0].position - myIntegrator.particleCompounds[0].position;
    relPosRef = msmrdtools::rotateVec(relPosRef, myIntegrator.particleCompounds[0].orientation.conj());
    REQUIRE(relPosRef == -1 * myIntegrator.pentamerCenter);

    // Check binding another particle into existing compound: bind particle 1 and 2 with bound state 1
    iIndex = 1;
    jIndex = 2;
    boundState = 1;
    myIntegrator.addCompound(plist, iIndex, jIndex, boundState);
    REQUIRE(myIntegrator.particleCompounds.size() == 1);
    REQUIRE(myIntegrator.getCompoundSize(0) == 3);
    REQUIRE(myIntegrator.particleCompounds[0].relativePositions.size() == 3);
    REQUIRE(myIntegrator.particleCompounds[0].relativeOrientations.size() == 3);
    // Check particles are set into the correct particleCompound index
    REQUIRE(plist[0].compoundIndex == 0);
    REQUIRE(plist[1].compoundIndex == 0);
    REQUIRE(plist[2].compoundIndex == 0);
    relPosition = myIntegrator.discreteTrajClass->getRelativePosition(boundState);
    relOrientation = myIntegrator.discreteTrajClass->getRelativeOrientation(boundState);
    newRelPosition = msmrdtools::rotateVec(relPosition, plist[iIndex].orientation);
    REQUIRE((plist[iIndex].position + newRelPosition - plist[jIndex].position).norm() <= 0.000001);
    REQUIRE((plist[iIndex].orientation * relOrientation - -1*plist[jIndex].orientation).norm() <= 0.000001);
    // Check relative position and orientation of previous binding is conserved
    iIndex = 0;
    jIndex = 2;
    boundState = 0;
    relPosition = myIntegrator.discreteTrajClass->getRelativePosition(boundState);
    relOrientation = myIntegrator.discreteTrajClass->getRelativeOrientation(boundState);
    newRelPosition = msmrdtools::rotateVec(relPosition, plist[iIndex].orientation);
    REQUIRE(plist[iIndex].position + newRelPosition == plist[jIndex].position);
    REQUIRE(plist[iIndex].orientation * relOrientation == plist[jIndex].orientation);

    // Check relative positions are stored correctly
    double angleDiff = 3 * M_PI / 5;
    vec3<double> relPos1 = {std::cos(angleDiff / 2.0), std::sin(angleDiff / 2.0), 0};
    vec3<double> relPos2 = {std::cos(angleDiff / 2.0), std::sin(-angleDiff / 2.0), 0};
    vec3<double> relPos1orthogonal = {-1.0 * std::sin(angleDiff / 2.0), std::cos(angleDiff / 2.0), 0.0};
    vec3<double> relPos2orthogonal = {std::sin(angleDiff / 2.0), std::cos(angleDiff / 2.0), 0.0};
    std::array<vec3<double>, 4> rotations;
    rotations[0] = {0.0, 0.0, -2 * M_PI / 5.0}; // part1Patch1 with part2patch1
    rotations[1] = M_PI * relPos1orthogonal; // part1Patch1 with part2patch2
    // --first 2 rotations correspond to binding on top patch (1) of part1, next 2 rotations to bottom patch (2).
    rotations[2] = {0.0, 0.0, 2 * M_PI / 5.0}; // part1Patch2 with part2patch1
    rotations[3] = M_PI * relPos2orthogonal; // part1Patch2 with part2patch2
    /*Convert rotations in the axis angle representation to quaternions */
    std::array<quaternion<double>, 4> quatRotations;
    for (int i = 0; i < 4; i++) {
        quatRotations[i] = msmrdtools::axisangle2quaternion(rotations[i]);
    }
    REQUIRE(-1*myIntegrator.pentamerCenter == myIntegrator.particleCompounds[0].relativePositions[0]);
    REQUIRE(-1*myIntegrator.pentamerCenter + relPos1 == myIntegrator.particleCompounds[0].relativePositions[2]);
    REQUIRE(-1*myIntegrator.pentamerCenter + relPos1 + msmrdtools::rotateVec(relPos1, quatRotations[0]) ==
             myIntegrator.particleCompounds[0].relativePositions[1]);

    // Make another independent complex: bind particle 3 and 4 with bound state 2
    iIndex = 3;
    jIndex = 4;
    boundState = 2;
    myIntegrator.addCompound(plist, iIndex, jIndex, boundState);
    REQUIRE(myIntegrator.particleCompounds.size() == 2);
    REQUIRE(myIntegrator.getCompoundSize(1) == 2);
    // Check particles are set into the correct particleCompound index
    REQUIRE(plist[iIndex].compoundIndex == 1);
    REQUIRE(plist[jIndex].compoundIndex == 1);
    // Check the newly bound particle is set into the correct relative position and orientation
    relPosition = myIntegrator.discreteTrajClass->getRelativePosition(boundState);
    relOrientation = myIntegrator.discreteTrajClass->getRelativeOrientation(boundState);
    newRelPosition = msmrdtools::rotateVec(relPosition, plist[iIndex].orientation);
    REQUIRE(plist[iIndex].position + newRelPosition == plist[jIndex].position);
    REQUIRE(plist[iIndex].orientation * relOrientation == plist[jIndex].orientation);

    // Diffuse compounds
    int timesteps = 20;
    auto prevPosition = 1.0 * myIntegrator.particleCompounds[0].position;
    auto prevOrientation = 1.0 * myIntegrator.particleCompounds[0].orientation;
    for (int i = 0; i < timesteps; i++) {
        myIntegrator.integrateDiffusionCompounds(plist,0.0001);
    }
    REQUIRE(prevPosition != myIntegrator.particleCompounds[0].position);
    REQUIRE(prevOrientation != myIntegrator.particleCompounds[0].orientation);

    // Check relative position and orientation is maintained after integration
    newRelPosition = msmrdtools::rotateVec(relPosition, plist[iIndex].orientation);
    REQUIRE((plist[iIndex].position + newRelPosition - plist[jIndex].position).norm() <= 0.000001);
    REQUIRE((plist[iIndex].orientation * relOrientation - plist[jIndex].orientation).norm() <= 0.000001);
    // Also check for other compound with two bonds (0-2 and 2-1)
    relPosition = myIntegrator.discreteTrajClass->getRelativePosition(0);
    relOrientation = myIntegrator.discreteTrajClass->getRelativeOrientation(0);
    newRelPosition = msmrdtools::rotateVec(relPosition, plist[0].orientation);
    REQUIRE((plist[0].position + newRelPosition - plist[2].position).norm() <= 0.000001);
    REQUIRE((plist[0].orientation * relOrientation - plist[2].orientation).norm() <= 0.000001);
    relPosition = myIntegrator.discreteTrajClass->getRelativePosition(1);
    relOrientation = myIntegrator.discreteTrajClass->getRelativeOrientation(1);
    newRelPosition = msmrdtools::rotateVec(relPosition, plist[2].orientation); // THIS IS WRONG BUT WHY??;
    REQUIRE((plist[2].position + newRelPosition - plist[1].position).norm() <= 0.000001);
    REQUIRE((plist[2].orientation * relOrientation - plist[1].orientation).norm() <= 0.000001);

    // Bind two existing compounds together (compound 0-2-1 with compound 3-4)
    iIndex = 0;
    jIndex = 3;
    boundState = 2;
    myIntegrator.addCompound(plist, iIndex, jIndex, boundState);
    REQUIRE(myIntegrator.particleCompounds.size() == 2);
    REQUIRE(myIntegrator.getCompoundSize(0) == 5);
    REQUIRE(myIntegrator.getCompoundSize(1) == 0);
    myIntegrator.cleanParticleCompoundsVector(plist);
    REQUIRE(myIntegrator.particleCompounds.size() == 1);
    REQUIRE(myIntegrator.getCompoundSize(0) == 5);
    //REQUIRE(myIntegrator.getCompoundSize(0) == 5);
    // Check relative positions/orientations match pentamer ring.
    // New binding 0-3 with bound state 3
    relPosition = myIntegrator.discreteTrajClass->getRelativePosition(boundState);
    relOrientation = myIntegrator.discreteTrajClass->getRelativeOrientation(boundState);
    newRelPosition = msmrdtools::rotateVec(relPosition, plist[iIndex].orientation);
    REQUIRE(plist[iIndex].position + newRelPosition == plist[jIndex].position);
    REQUIRE(plist[iIndex].orientation * relOrientation == plist[jIndex].orientation);
    // Check binding 0-2 with bound state 0
    relPosition = myIntegrator.discreteTrajClass->getRelativePosition(0);
    relOrientation = myIntegrator.discreteTrajClass->getRelativeOrientation(0);
    newRelPosition = msmrdtools::rotateVec(relPosition, plist[0].orientation);
    REQUIRE((plist[0].position + newRelPosition - plist[2].position).norm() <= 0.000001);
    REQUIRE((plist[0].orientation * relOrientation- plist[2].orientation).norm() <= 0.000001);
    // Check binding 2-1 with bound state 1
    relPosition = myIntegrator.discreteTrajClass->getRelativePosition(1);
    relOrientation = myIntegrator.discreteTrajClass->getRelativeOrientation(1);
    newRelPosition = msmrdtools::rotateVec(relPosition, plist[2].orientation);
    REQUIRE((plist[2].position + newRelPosition - plist[1].position).norm() <= 0.000001);
    REQUIRE((plist[2].orientation * relOrientation - plist[1].orientation).norm() <= 0.000001);
    // Check binding 3-4 with bound state 2
    relPosition = myIntegrator.discreteTrajClass->getRelativePosition(2);
    relOrientation = myIntegrator.discreteTrajClass->getRelativeOrientation(2);
    newRelPosition = msmrdtools::rotateVec(relPosition, plist[3].orientation);
    REQUIRE((plist[3].position + newRelPosition - plist[4].position).norm() <= 0.000001);
    REQUIRE((plist[3].orientation * relOrientation - plist[4].orientation).norm() <= 0.000001);

    // Close compound into a pentamer (closes it by binding 4 with 1) Closing compound deactivated for now.
    iIndex = 1;
    jIndex = 4;
    boundState = 2;
    myIntegrator.addCompound(plist, iIndex, jIndex, boundState);
    REQUIRE(myIntegrator.getNumberOfBindingsInCompound(0) == 5);
    auto bindingLoops = myIntegrator.findClosedBindingLoops();
    REQUIRE(bindingLoops[0] == 5);
}
