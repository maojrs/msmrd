//
// Created by maojrs on 3/28/19.
//
#include <catch2/catch.hpp>
#include "discretizations/spherePartition.hpp"
#include "tools.hpp"

using namespace msmrd;

TEST_CASE("Spherical partition", "[spherePartition]") {
int numSections = 15;
auto spherePart = spherePartition(numSections);
auto regionsPerCollar = spherePart.regionsPerCollar;
auto phis = spherePart.phis;
auto thetas = spherePart.thetas;
std::vector<int> regionsPerCollarRef{1, 6, 7, 1};
std::vector<double> phisRef{0.0, 0.52231482, 1.50408018, 2.61927783};
std::vector<std::vector<double>> thetasRef;
thetasRef.resize(2);
thetasRef[0] = std::vector<double> {0.0, 1.04719755, 2.0943951 , 3.14159265, 4.1887902 , 5.23598776};
thetasRef[1] = std::vector<double> {0.0, 0.8975979 , 1.7951958 , 2.6927937 , 3.5903916 , 4.48798951, 5.38558741};
REQUIRE(regionsPerCollar == regionsPerCollarRef);
REQUIRE(msmrdtools::stdvecNorm(phis, phisRef) <= 0.000001);
REQUIRE(msmrdtools::stdvecNorm(thetas[0], thetasRef[0]) <= 0.000001);
REQUIRE(msmrdtools::stdvecNorm(thetas[1], thetasRef[1]) <= 0.000001);
// Now test getSectionNumber function
vec3<double> coordinateUp{0.0, 0.0, 1.0};
vec3<double> coordinateDown{0.0, 0.0, -1.0};
vec3<double> coordinate1{1.0, 2.0, 0.0};
vec3<double> coordinate2{0.5, 3.5, 0.3};
vec3<double> coordinate3{-2.5, -1, -0.8};
int secNumUp = spherePart.getSectionNumber(coordinateUp);
int secNumDown = spherePart.getSectionNumber(coordinateDown);
int secNum1 = spherePart.getSectionNumber(coordinate1);
int secNum2 = spherePart.getSectionNumber(coordinate2);
int secNum3 = spherePart.getSectionNumber(coordinate3);
REQUIRE(secNumUp == 1);
REQUIRE(secNumDown == 15);
REQUIRE(secNum1 == 9);
REQUIRE(secNum2 == 3);
REQUIRE(secNum3 == 11);
// Now test getAngles function
auto anglesUp = spherePart.getAngles(secNumUp);
auto anglesDown = spherePart.getAngles(secNumDown);
auto angles1 = spherePart.getAngles(secNum1);
auto angles2 = spherePart.getAngles(secNum2);
auto angles3 = spherePart.getAngles(secNum3);
auto phiIntervalUp = std::get<0>(anglesUp);
auto thetaIntervalUp = std::get<1>(anglesUp);
auto phiIntervalDown = std::get<0>(anglesDown);
auto thetaIntervalDown = std::get<1>(anglesDown);
auto phiInterval1 = std::get<0>(angles1);
auto thetaInterval1 = std::get<1>(angles1);
auto phiInterval2 = std::get<0>(angles2);
auto thetaInterval2 = std::get<1>(angles2);
auto phiInterval3 = std::get<0>(angles3);
auto thetaInterval3 = std::get<1>(angles3);
REQUIRE(phiIntervalUp == std::vector<double>{0.0, 0.5223148218060486} );
REQUIRE(thetaIntervalUp == std::vector<double>{0, 6.283185307179586} );
REQUIRE(phiIntervalDown == std::vector<double>{2.6192778317837444, 3.141592653589793} );
REQUIRE(thetaIntervalDown == std::vector<double>{0, 6.283185307179586} );
REQUIRE(phiInterval1 == std::vector<double>{1.5040801783846711, 2.6192778317837444} );
REQUIRE(thetaInterval1 == std::vector<double>{0.8975979010256552, 1.7951958020513104} );
REQUIRE(phiInterval2 == std::vector<double>{0.5223148218060486, 1.5040801783846711} );
REQUIRE(thetaInterval2 == std::vector<double>{1.0471975511965976, 2.0943951023931953} );
REQUIRE(phiInterval3 == std::vector<double>{1.5040801783846711, 2.6192778317837444} );
REQUIRE(thetaInterval3 == std::vector<double>{2.6927937030769655, 3.5903916041026207} );
}