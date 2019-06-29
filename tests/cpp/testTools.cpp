//
// Created by maojrs on 5/29/19.
//
#include <catch2/catch.hpp>
#include <randomgen.hpp>
#include "particle.hpp"
#include "tools.hpp"
#include "vec3.hpp"

using namespace msmrd;


TEST_CASE("Relative distance calculation", "[calculateRelativePosition]") {
    vec3<double> v1({-0.2, -0.3, 0.4});
    vec3<double> v2({0.4, 0.3, -0.1});
    bool boundaryActive = true;
    std::string boundaryType = "periodic";
    vec3<double> boxsize = {1., 1., 1.};
    vec3<double> relPosition;
    relPosition = msmrdtools::calculateRelativePosition(v1, v2, boundaryActive, boundaryType, boxsize);
    vec3<double> refSolution = {-0.4, -0.4, -0.5};
    REQUIRE( (relPosition -refSolution).norm() < 0.0000001 );
}

TEST_CASE("Switching between quaternion and axis angle representations", "[quaternion2axisangle, "
                                                                         "axisangle2quaternion]") {
    randomgen randg = randomgen();
    double theta = randg.uniformRange(0, 2*M_PI);
    vec3<double> axisAngleRep = randg.normal3D(0,5);
    axisAngleRep = theta*axisAngleRep/axisAngleRep.norm();
    quaternion<double> quatRep = msmrdtools::axisangle2quaternion(axisAngleRep);
    vec3<double> recoveredAxisAngleRep = msmrdtools::quaternion2axisangle(quatRep);
    REQUIRE((axisAngleRep - recoveredAxisAngleRep).norm() < 0.000001);
}
