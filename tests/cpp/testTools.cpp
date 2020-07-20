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

TEST_CASE("Rotate vector through an axis off the origin", "[rotateVecOffAxis]") {
    vec3<double> vec1 = vec3<double>(1,1,0);
    vec3<double> offAxisPoint = vec3<double>(2,0,0);
    vec3<double> rotation = vec3<double>(0,0,M_PI);
    auto quatRotation = msmrdtools::axisangle2quaternion(rotation);
    auto resultVec = msmrdtools::rotateVecOffAxis(vec1, quatRotation,offAxisPoint);
    auto referenceResult = vec3<double>(3,-1,0);
    REQUIRE((resultVec - referenceResult).norm() < 0.000001);
}

TEST_CASE("Recover rotation quaternion from vectors", "[recoverRotationFromVectors]") {
    // Simple but tricky rotation (rotation by pi yields zero cross product)
    vec3<double> origin = vec3<double>(0,0,0);
    vec3<double> vec1 = vec3<double>(1,1,0);
    vec3<double> vec2 = vec3<double>(1,-1,0);
    vec3<double> offAxisPoint = vec3<double>(2,0,0);
    vec3<double> rotation = vec3<double>(0,0,M_PI);
    auto quatRotation = msmrdtools::axisangle2quaternion(rotation);
    auto newOrigin = msmrdtools::rotateVecOffAxis(origin, quatRotation, offAxisPoint);
    auto rotatedVec1 = msmrdtools::rotateVecOffAxis(vec1, quatRotation, offAxisPoint);
    auto rotatedVec2 = msmrdtools::rotateVecOffAxis(vec2, quatRotation, offAxisPoint);
    auto recoveredQuaternion = msmrdtools::recoverRotationFromVectors(origin,vec1,vec2,
            newOrigin,rotatedVec1,rotatedVec2);
    auto recoveredAxisAngle = msmrdtools::quaternion2axisangle(recoveredQuaternion);
    REQUIRE((rotation - recoveredAxisAngle).norm() < 0.000001);
    REQUIRE((quatRotation - recoveredQuaternion).norm() < 0.000001);
    // Very arbitrary rotation
    rotation = vec3<double>(0.4,0.7*M_PI,0.2*M_PI);
    quatRotation = msmrdtools::axisangle2quaternion(rotation);
    newOrigin = msmrdtools::rotateVecOffAxis(origin, quatRotation, offAxisPoint);
    rotatedVec1 = msmrdtools::rotateVecOffAxis(vec1, quatRotation, offAxisPoint);
    rotatedVec2 = msmrdtools::rotateVecOffAxis(vec2, quatRotation, offAxisPoint);
    recoveredQuaternion = msmrdtools::recoverRotationFromVectors(origin,vec1,vec2,
            newOrigin,rotatedVec1,rotatedVec2);
    recoveredAxisAngle = msmrdtools::quaternion2axisangle(recoveredQuaternion);
    REQUIRE((rotation - recoveredAxisAngle).norm() < 0.000001);
    REQUIRE((quatRotation - recoveredQuaternion).norm() < 0.000001);


}
