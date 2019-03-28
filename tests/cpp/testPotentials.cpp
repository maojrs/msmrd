//
// Created by maojrs on 3/28/19.
//

#include <catch2/catch.hpp>
#include "potentials/gayBerne.hpp"
#include "potentials/patchyParticle.hpp"

using namespace msmrd;

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
