//
// Created by maojrs on 3/28/19.
//
#include <catch2/catch.hpp>
#include "discretizations/spherePartition.hpp"
#include "discretizations/halfSpherePartition.hpp"
#include "discretizations/quaternionPartition.hpp"
#include "discretizations/positionOrientationPartition.hpp"
#include "tools.hpp"

using namespace msmrd;


TEST_CASE("Spherical partition", "[spherePartition]") {
    // Test main partition function
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
    REQUIRE(phiIntervalUp == std::array<double,2>{0.0, 0.5223148218060486} );
    REQUIRE(thetaIntervalUp == std::array<double,2>{0, 6.283185307179586} );
    REQUIRE(phiIntervalDown == std::array<double,2>{2.6192778317837444, 3.141592653589793} );
    REQUIRE(thetaIntervalDown == std::array<double,2>{0, 6.283185307179586} );
    REQUIRE(phiInterval1 == std::array<double,2>{1.5040801783846711, 2.6192778317837444} );
    REQUIRE(thetaInterval1 == std::array<double,2>{0.8975979010256552, 1.7951958020513104} );
    REQUIRE(phiInterval2 == std::array<double,2>{0.5223148218060486, 1.5040801783846711} );
    REQUIRE(thetaInterval2 == std::array<double,2>{1.0471975511965976, 2.0943951023931953} );
    REQUIRE(phiInterval3 == std::array<double,2>{1.5040801783846711, 2.6192778317837444} );
    REQUIRE(thetaInterval3 == std::array<double,2>{2.6927937030769655, 3.5903916041026207} );
}


TEST_CASE("Half spherical partition", "[halfSpherePartition]") {
    // Test main partition function
    int numSections = 15;
    auto halfSpherePart = halfSpherePartition(numSections);
    auto regionsPerCollar = halfSpherePart.regionsPerCollar;
    auto phis = halfSpherePart.phis;
    auto thetas = halfSpherePart.thetas;
    std::vector<int> regionsPerCollarRef{1, 4, 5, 4, 1};
    std::vector<double> phisRef{0.0, 0.52231482, 1.23095942, 1.91063324, 2.61927783};
    std::vector<std::vector<double>> thetasRef;
    thetasRef.resize(3);
    thetasRef[0] = std::vector<double> {0.0, 0.78539816, 1.57079633, 2.35619449};
    thetasRef[1] = std::vector<double> {0.0, 0.62831853, 1.25663706, 1.88495559, 2.51327412};
    thetasRef[2] = std::vector<double> {0.0, 0.78539816, 1.57079633, 2.35619449};
    REQUIRE(regionsPerCollar == regionsPerCollarRef);
    REQUIRE(msmrdtools::stdvecNorm(phis, phisRef) <= 0.000001);
    REQUIRE(msmrdtools::stdvecNorm(thetas[0], thetasRef[0]) <= 0.000001);
    REQUIRE(msmrdtools::stdvecNorm(thetas[1], thetasRef[1]) <= 0.000001);
    REQUIRE(msmrdtools::stdvecNorm(thetas[2], thetasRef[2]) <= 0.000001);

    // Now test getSectionNumber function
    vec3<double> coordinateUp{0.0, 0.0, 1.0};
    vec3<double> coordinateDown{0.0, 0.0, -1.0};
    vec3<double> coordinate1{0.0, 1.0, 0.0};
    vec3<double> coordinate2{0.9, 1.0, 1.0};
    vec3<double> coordinate3{-0.9, 1, -1.0};
    int secNumUp = halfSpherePart.getSectionNumber(coordinateUp);
    int secNumDown = halfSpherePart.getSectionNumber(coordinateDown);
    int secNum1 = halfSpherePart.getSectionNumber(coordinate1);
    int secNum2 = halfSpherePart.getSectionNumber(coordinate2);
    int secNum3 = halfSpherePart.getSectionNumber(coordinate3);
    REQUIRE(secNumUp == 1);
    REQUIRE(secNumDown == 15);
    REQUIRE(secNum1 == 8);
    REQUIRE(secNum2 == 3);
    REQUIRE(secNum3 == 13);

    // Now test getAngles function
    auto anglesUp = halfSpherePart.getAngles(1);
    auto anglesDown = halfSpherePart.getAngles(numSections);
    auto angles1 = halfSpherePart.getAngles(4);
    auto angles2 = halfSpherePart.getAngles(7);
    auto angles3 = halfSpherePart.getAngles(11);
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
    REQUIRE(phiIntervalUp == std::array<double,2>{0.0, 0.5223148218060486} );
    REQUIRE(thetaIntervalUp == std::array<double,2>{0, 3.141592653589793} );
    REQUIRE(phiIntervalDown == std::array<double,2>{2.6192778317837444, 3.141592653589793} );
    REQUIRE(thetaIntervalDown == std::array<double,2>{0, 3.141592653589793} );
    REQUIRE(phiInterval1 == std::array<double,2>{0.5223148218060486, 1.2309594173407745} );
    REQUIRE(thetaInterval1 == std::array<double,2>{1.5707963267948966, 2.356194490192345} );
    REQUIRE(phiInterval2 == std::array<double,2>{1.2309594173407745, 1.9106332362490186} );
    REQUIRE(thetaInterval2 == std::array<double,2>{0.6283185307179586, 1.2566370614359172} );
    REQUIRE(phiInterval3 == std::array<double,2>{1.9106332362490186, 2.6192778317837444} );
    REQUIRE(thetaInterval3 == std::array<double,2>{0.0, 0.7853981633974483} );
}


TEST_CASE("Quaternion partition", "[quaternionPartition]") {
    // Test main partition function
    int numRadialSections = 5;
    int numSphericalSections = 15;
    auto quatPartition = quaternionPartition(numRadialSections, numSphericalSections);
    auto rslices = quatPartition.radialSections;
    auto regionsPerCollar = quatPartition.sphericalPartition->regionsPerCollar;
    auto phis = quatPartition.sphericalPartition->phis;
    auto thetas = quatPartition.sphericalPartition->thetas;
    std::vector<double> rslicesRef{ 0.0, 0.2, 0.4, 0.6, 0.8, 1.0 };
    std::vector<int> regionsPerCollarRef{1, 6, 7, 1};
    std::vector<double> phisRef{0.0, 0.52231482, 1.50408018, 2.61927783};
    std::vector<std::vector<double>> thetasRef;
    thetasRef.resize(2);
    thetasRef[0] = std::vector<double> {0.0, 1.04719755, 2.0943951 , 3.14159265, 4.1887902 , 5.23598776};
    thetasRef[1] = std::vector<double> {0.0, 0.8975979 , 1.7951958 , 2.6927937 , 3.5903916 , 4.48798951, 5.38558741};
    REQUIRE(msmrdtools::stdvecNorm(rslices, rslicesRef) <= 0.000001);
    REQUIRE(regionsPerCollar == regionsPerCollarRef);
    REQUIRE(msmrdtools::stdvecNorm(phis, phisRef) <= 0.000001);
    REQUIRE(msmrdtools::stdvecNorm(thetas[0], thetasRef[0]) <= 0.000001);
    REQUIRE(msmrdtools::stdvecNorm(thetas[1], thetasRef[1]) <= 0.000001);

    // Test getSectionNumber from coordinate (s,x,z)
    quaternion<double> coordinate1{0.0, 0.0, 0.0, 0.1};
    quaternion<double> coordinate2{0.0, 0.0, 0.0, 0.25};
    quaternion<double> coordinate3{0.0, 0.0, 0.0, -0.25};
    quaternion<double> coordinate4{0.0, 0.0, 0.0, 0.85};
    quaternion<double> coordinate5{0.0, 0.0, 0.0, -0.85};
    quaternion<double> coordinate6{0.0, 1.0, 2.0, 0.0};
    quaternion<double> coordinate7{-1.0, -0.5, -3.5, -0.3};
    quaternion<double> coordinate8{0.2, -2.5, -1, -0.8};
    coordinate6 = 0.65*coordinate6/coordinate6.norm();
    coordinate7 = 0.43*coordinate7/coordinate7.norm();
    coordinate8 = 0.82*coordinate8/coordinate8.norm();
    int secNum1 = quatPartition.getSectionNumber(coordinate1);
    int secNum2 = quatPartition.getSectionNumber(coordinate2);
    int secNum3 = quatPartition.getSectionNumber(coordinate3);
    int secNum4 = quatPartition.getSectionNumber(coordinate4);
    int secNum5 = quatPartition.getSectionNumber(coordinate5);
    int secNum6 = quatPartition.getSectionNumber(coordinate6);
    int secNum7 = quatPartition.getSectionNumber(coordinate7);
    int secNum8 = quatPartition.getSectionNumber(coordinate8);
    REQUIRE(secNum1 == 1);
    REQUIRE(secNum2 == 2);
    REQUIRE(secNum3 == 16);
    REQUIRE(secNum4 == 47);
    REQUIRE(secNum5 == 61);
    REQUIRE(secNum6 == 31 + 9);
    REQUIRE(secNum7 == 16 + 3);
    REQUIRE(secNum8 == 46 + 11);

    // Test getSectionIntervals from sectionNumber
    auto intervals0 = quatPartition.getSectionIntervals(1); // Special (whole inner sphere)
    auto intervals1 = quatPartition.getSectionIntervals(61);
    auto intervals2 = quatPartition.getSectionIntervals(40);
    auto intervals3 = quatPartition.getSectionIntervals(19);
    auto intervals4 = quatPartition.getSectionIntervals(57);
    // r intervals
    std::array<double, 2> rIntervalRef0 = {0.0, 0.2};
    std::array<double, 2> rIntervalRef1 = {0.8, 1.0};
    std::array<double, 2> rIntervalRef2 = {0.6, 0.8};
    std::array<double, 2> rIntervalRef3 = {0.4, 0.6};
    std::array<double, 2> rIntervalRef4 = {0.8, 1.0};
    REQUIRE(msmrdtools::stdvecNorm(std::get<0>(intervals0), rIntervalRef0) <= 0.000001);
    REQUIRE(msmrdtools::stdvecNorm(std::get<0>(intervals1), rIntervalRef1) <= 0.000001);
    REQUIRE(msmrdtools::stdvecNorm(std::get<0>(intervals2), rIntervalRef2) <= 0.000001);
    REQUIRE(msmrdtools::stdvecNorm(std::get<0>(intervals3), rIntervalRef3) <= 0.000001);
    REQUIRE(msmrdtools::stdvecNorm(std::get<0>(intervals4), rIntervalRef4) <= 0.000001);
    // angular intervals (the -1 correspond to the inner sphere)
    auto anglesRef1 = quatPartition.sphericalPartition->getAngles(15); // (61-1)%15 = 0 -> 15
    auto anglesRef2 = quatPartition.sphericalPartition->getAngles((40-1)%15);
    auto anglesRef3 = quatPartition.sphericalPartition->getAngles((19-1)%15);
    auto anglesRef4 = quatPartition.sphericalPartition->getAngles((57-1)%15);
    // phiIntervals
    std::array<double, 2> phiRef0 = {0, M_PI};
    REQUIRE(msmrdtools::stdvecNorm(std::get<1>(intervals0), phiRef0) <= 0.000001);
    REQUIRE(msmrdtools::stdvecNorm(std::get<1>(intervals1), std::get<0>(anglesRef1)) <= 0.000001);
    REQUIRE(msmrdtools::stdvecNorm(std::get<1>(intervals2), std::get<0>(anglesRef2)) <= 0.000001);
    REQUIRE(msmrdtools::stdvecNorm(std::get<1>(intervals3), std::get<0>(anglesRef3)) <= 0.000001);
    REQUIRE(msmrdtools::stdvecNorm(std::get<1>(intervals4), std::get<0>(anglesRef4)) <= 0.000001);
    // theta intervals
    std::array<double, 2> thetaRef0 = {0, 2*M_PI};
    REQUIRE(msmrdtools::stdvecNorm(std::get<2>(intervals0), thetaRef0) <= 0.000001);
    REQUIRE(msmrdtools::stdvecNorm(std::get<2>(intervals1), std::get<1>(anglesRef1)) <= 0.000001);
    REQUIRE(msmrdtools::stdvecNorm(std::get<2>(intervals2), std::get<1>(anglesRef2)) <= 0.000001);
    REQUIRE(msmrdtools::stdvecNorm(std::get<2>(intervals3), std::get<1>(anglesRef3)) <= 0.000001);
    REQUIRE(msmrdtools::stdvecNorm(std::get<2>(intervals4), std::get<1>(anglesRef4)) <= 0.000001);
}

TEST_CASE("position orientation partition", "[positionOrientationPartition]") {
    // Create position orientation partition (six dimensional)
    int numSphericalSectionsPos = 7;
    int numRadialSectionsQuat = 5;
    int numSphericalSectionsQuat = 7;
    auto positionOrientationPart = new positionOrientationPartition(2.2,
            numSphericalSectionsPos, numRadialSectionsQuat, numSphericalSectionsQuat);

    int totalnumSecsQuat = numSphericalSectionsQuat*(numRadialSectionsQuat -1) + 1;
    int numTransitionsStates = numSphericalSectionsPos * totalnumSecsQuat; //203

    std::vector<int> endStates = {0, 50, 100, 150, numTransitionsStates};
    for (auto const &endState : endStates) {
        auto intervals = positionOrientationPart->getSectionIntervals(endState);
        auto interval = std::get<0>(intervals);
        REQUIRE(0 <= interval[0]);
        REQUIRE(interval[0] < interval[1]);
        interval = std::get<1>(intervals);
        REQUIRE(0 <= interval[0]);
        REQUIRE(interval[0] <= interval[1]);
        interval = std::get<2>(intervals);
        REQUIRE(0 <= interval[0]);
        REQUIRE(interval[0] <= interval[1]);
        interval = std::get<3>(intervals);
        REQUIRE(0 <= interval[0]);
        REQUIRE(interval[0] <= interval[1]);
        interval = std::get<4>(intervals);
        REQUIRE(0 <= interval[0]);
        REQUIRE(interval[0] <= interval[1]);
    }
    auto numTotalSecs = positionOrientationPart->numTotalSections;
    REQUIRE(numTotalSecs == 203);
}
