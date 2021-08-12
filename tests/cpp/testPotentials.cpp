//
// Created by maojrs on 3/28/19.
//

#include <catch2/catch.hpp>
#include "randomgen.hpp"
#include "potentials/gaussians3D.hpp"
#include "potentials/gayBerne.hpp"
#include "potentials/patchyParticle.hpp"
#include "potentials/patchyProteinMarkovSwitch.hpp"

using namespace msmrd;

TEST_CASE("External Potentials (Gaussians3D))", "[potentials]") {
    std::array<vec3<double>, 2> forctorq1;
    std::array<vec3<double>, 2> forctorq2;

    // Create Gaussian potential with one minima
    std::vector<std::vector<double>> minimaPositions;
    std::vector<std::vector<double>> stdDeviations;
    std::vector<int> partTypes;
    double scalefactor = 500;

    minimaPositions.resize(1);
    stdDeviations.resize(1);
    minimaPositions[0] = std::vector<double>{0,0,0};
    stdDeviations[0] = std::vector<double>{5,5,5};
    partTypes = std::vector<int>{1};
    auto potentialGauss = gaussians3D(minimaPositions, stdDeviations, partTypes, scalefactor);
    REQUIRE(potentialGauss.nminima == 1);

    // Create particle to feel potential (one is type one other one is not, potential should only act on type one)
    vec3<double> pos = vec3<double>(3.,2.,1.);
    vec3<double> theta = vec3<double>(0., 1., 0);
    particle part1 = particle(1.0, 1.0, pos, theta);
    particle part2 = particle(1.0, 1.0, pos, theta);
    part1.setType(1);

    // Test potentials
    forctorq1 = potentialGauss.forceTorque(part1);
    forctorq2 = potentialGauss.forceTorque(part2);
    REQUIRE(forctorq1[0] != vec3<double>{0,0,0});
    REQUIRE(forctorq2[0] == vec3<double>{0,0,0});
}


TEST_CASE("Pair Potentials consistency", "[potentials]") {
    std::array<vec3<double>, 4> forctorq1;
    std::array<vec3<double>, 4> forctorq2;

    // Gay Berne potential consistency check
    auto potentialGB = gayBerne(3.0,5.0,5.0,1.0);
    double th = 1.1*M_PI/2.0;
    vec3<double> pos1 = vec3<double>(0.,0.,0.);
    vec3<double> pos2 = vec3<double>(1.,0.,0.);
    vec3<double> th1 = vec3<double>(0., 1., 0);
    vec3<double> th2 = vec3<double>(std::cos(th), std::sin(th),0.);
    particle part1 = particle(1.0, 1.0, pos1, th1);
    particle part2 = particle(1.0, 1.0, pos2, th2);
    forctorq1 = potentialGB.forceTorque(part1, part2);
    forctorq2 = potentialGB.forceTorque(part2, part1);
    // Check consistency between forces
    REQUIRE(forctorq1[0] == forctorq2[2]);
    REQUIRE(forctorq2[0] == forctorq1[2]);
    // Check consistency between torques
    REQUIRE(forctorq1[1] == forctorq2[3]);
    REQUIRE(forctorq2[1] == forctorq1[3]);

    // Patchy particle potential consistency check
    double sigma = 1.0;
    double strength = 100.0;
    auto patch1 = vec3<double> (1.,0.,0.);
    std::vector<vec3<double>> patchesCoordinates = {patch1};
    auto potentialPatchyParticle = patchyParticle(sigma, strength, patchesCoordinates);
    th = 0.9*M_PI;
    pos1 = vec3<double>(0.,0.,0.);
    pos2 = vec3<double>(1.0,-0.3,0.);
    quaternion<double> q1 = quaternion<double>(1.,0.,0., 0.);
    quaternion<double> q2 = quaternion<double>(std::cos(th/2.0), 0., 0., std::sin(th/2.0));
    part1 = particle(1.0, 1.0, pos1, q1);
    part2 = particle(1.0, 1.0, pos2, q2);
    forctorq1 = potentialPatchyParticle.forceTorque(part1, part2);
    forctorq2 = potentialPatchyParticle.forceTorque(part2, part1);
    // Check consistency between forces
    REQUIRE(forctorq1[0] == forctorq2[2]);
    REQUIRE(forctorq2[0] == forctorq1[2]);
    // Check consistency between torques
    REQUIRE(forctorq1[1] == forctorq2[3]);
    REQUIRE(forctorq2[1] == forctorq1[3]);
}

TEST_CASE("patchyProteinMS potential: test calculatePlanes function", "[potentials]") {
    // Define patchy protein potential
    std::vector<vec3<double>> patchesCoordinatesA(6);
    std::vector<vec3<double>> patchesCoordinatesB(1);
    patchesCoordinatesA[0] = vec3<double> (1.,0.,0.);
    patchesCoordinatesA[1] = vec3<double> (0.,1.,0.);
    patchesCoordinatesA[2] = vec3<double> (0.,0.,1.);
    patchesCoordinatesA[3] = vec3<double> (-1.,0.,0.);
    patchesCoordinatesA[4] = vec3<double> (0.,-1.,0.);
    patchesCoordinatesA[5] = vec3<double> (0.,0.,-1.);
    patchesCoordinatesB[0] = vec3<double> (1.,0.,0.);
    double sigma = 1.0;
    double strength = 1.0;
    auto potentialPPMS = patchyProteinMarkovSwitch(sigma, strength, strength, patchesCoordinatesA, patchesCoordinatesB);
    // Define relative positions and orientations
    std::vector<vec3<double>> relativePositions(6);
    std::vector<vec3<double>> relativeOrientations(6);
    relativePositions[0] = {1., 0., 0.};
    relativePositions[1] = {0., 1., 0.};
    relativePositions[2] = {0., 0., 1.};
    relativePositions[3] = {-1., 0., 0.};
    relativePositions[4] = {0., -1., 0.};
    relativePositions[5] = {0., 0., -1.};
    relativeOrientations[0] = {0.0, 0.0, M_PI};
    relativeOrientations[1] = {0.0, 0.0, -M_PI / 2.0};
    relativeOrientations[2] = {0.0, M_PI / 2.0, 0.0};
    relativeOrientations[3] = {0.0, 0.0, 0.0};
    relativeOrientations[4] = {0.0, 0.0, M_PI / 2.0};
    relativeOrientations[5] = {0.0, -M_PI / 2.0, 0.0};
    // Define particles variables
    vec3<double> pos1;
    vec3<double> pos2;
    quaternion<double> th1;
    quaternion<double> th2;
    // Check output planes are the same if orientation is exact
    for (int i = 0; i<6; i++) {
        pos1 = vec3<double>(0.,0.,0.);
        pos2 = relativePositions[i];
        th1 = quaternion<double>(1., 0., 0., 0.);
        th2 = msmrdtools::axisangle2quaternion(relativeOrientations[i]);
        auto part1 = particle(0, 0, 1.0, 1.0, pos1, th1);
        auto part2 = particle(0, 0, 1.0, 1.0, pos2, th2);
        auto planes = potentialPPMS.calculatePlanes(part1, part2, patchesCoordinatesA, patchesCoordinatesB);
        auto plane1 = std::get<0>(planes);
        auto plane2 = std::get<1>(planes);
        REQUIRE((plane1 - plane2).norm() <= 1E-10);
        }
    // Check output planes are close if orientation is not exact
    randomgen randg = randomgen();
    for (int k = 0; k <1000; k++) {
        for (int i = 0; i < 6; i++) {
            pos1 = vec3<double>(0., 0., 0.) + randg.normal3D(0, 0.01);
            pos2 = relativePositions[i];
            th1 = quaternion<double>(1., 0., 0., 0.);
            auto relOrientations = relativeOrientations[i] + randg.normal3D(0, 0.01);
            th2 = msmrdtools::axisangle2quaternion(relOrientations);
            auto part1 = particle(0, 0, 1.0, 1.0, pos1, th1);
            auto part2 = particle(0, 0, 1.0, 1.0, pos2, th2);
            auto planes = potentialPPMS.calculatePlanes(part1, part2, patchesCoordinatesA, patchesCoordinatesB);
            auto plane1 = std::get<0>(planes);
            auto plane2 = std::get<1>(planes);
            REQUIRE((plane1 - plane2).norm() <= 0.1);
        }
    }
}

TEST_CASE("patchyProteinMS potential: test patchesActive", "[potentials]") {
    // Define patchy protein potential
    std::vector<vec3<double>> patchesCoordinatesA(6);
    std::vector<vec3<double>> patchesCoordinatesB(0);
    patchesCoordinatesA[0] = vec3<double> (1.,0.,0.);
    patchesCoordinatesA[1] = vec3<double> (0.,1.,0.);
    patchesCoordinatesA[2] = vec3<double> (0.,0.,1.);
    patchesCoordinatesA[3] = vec3<double> (-1.,0.,0.);
    patchesCoordinatesA[4] = vec3<double> (0.,-1.,0.);
    patchesCoordinatesA[5] = vec3<double> (0.,0.,-1.);
    double sigma = 1.0;
    double strength = 1.0;
    auto potentialPPMS = patchyProteinMarkovSwitch(sigma, strength, strength, patchesCoordinatesA, patchesCoordinatesB);
    REQUIRE(potentialPPMS.arePatchesActive() == false);
}

