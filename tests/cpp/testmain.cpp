//
// Created by dibakma on 26.06.18.
//
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include <catch2/catch.hpp>
#include "potentials/gayBerne.hpp"
#include "potentials/patchyParticle.hpp"
#include "markovModels/discreteTimeMarkovModel.hpp"
#include "markovModels/continuousTimeMarkovModel.hpp"
#include "integrators/overdampedLangevin.hpp"
#include "simulation.hpp"
#include "trajectories/trajectory.hpp"
#include "particle.hpp"
#include "quaternion.hpp"
#include "randomgen.hpp"
#include "spherePartition.hpp"
#include "tools.hpp"
#include "vec3.hpp"

using namespace msmrd;
using msm = msmrd::discreteTimeMarkovStateModel;
using ctmsm = msmrd::continuousTimeMarkovStateModel;

TEST_CASE("Basic vector arithmetic", "[vectors]") {
    vec3<double> v1({1.,2.,3});
    vec3<double> v2({3.,1.,2});
    vec3<double> v3({4.,3.,5});
    vec3<double> v4({2,4,6});
    REQUIRE( v1+v2 == v3 );
    REQUIRE( v1*v2 == 11);
    REQUIRE( 5*v1 == v1*5);
    REQUIRE( 2*v1 == v4 );
}

TEST_CASE("Basic quaternion arithmetic", "[quaternions]") {
    quaternion<double> q1({1.,2.,3.,4.});
    quaternion<double> q1p(1,2,3,4);
    quaternion<double> q2({5.,1.,7.,5.});
    quaternion<double> q3({6.,3.,10.,9.});
    std::array<double, 4> v1({5.,1.,7.,5.});
    std::array<double, 4> v3({6.,3.,10.,9.});
    std::array<double, 4> v4({-38,-2,16,36});
    REQUIRE( q1p == q1 );
    REQUIRE( q1+q2 == q3 );
    REQUIRE( q1+v1 == q3 );
    REQUIRE( q1p+v1 == v3 );
    REQUIRE( q1*q2 == v4 );
}

TEST_CASE("Sampling from randomgen", "[randomgen]") {
    randomgen randg;
    // Uniform Range
    double rand1;
    double rand2;
    double rand3;
    for (int i=0; i<100; i++) {
        rand1 = randg.uniformRange(-5,5);
        rand2 = randg.uniformRange(100,500);
        rand3 = randg.uniformRange(-750,-450);
        REQUIRE(rand1 >= -5);
        REQUIRE(rand1 <= 5);
        REQUIRE(rand2 >= 100);
        REQUIRE(rand2 <= 500);
        REQUIRE(rand3 >= -750);
        REQUIRE(rand3 <= -450);
    }
    // Normal distribution (might fail sporadically, maybe remove)
    rand1 = 0;
    for (int i=0; i<1000; i++) {
        rand1 += randg.normal(0,1)/1000;
    }
    REQUIRE(std::abs(rand1) <= 0.1);
    // Uniform Sphere
    vec3<double> trial1;
    vec3<double> trial2;
    vec3<double> trial3;
    for (int i=0; i<100; i++) {
        trial1 = randg.uniformSphere(1.0);
        trial2 = randg.uniformSphere(2.0);
        trial3 = randg.uniformSphere(3.0);
        REQUIRE(trial1.norm() <= 1.0);
        REQUIRE(trial2.norm() <= 2.0);
        REQUIRE(trial3.norm() <= 3.0);
    }
    // Uniform Shell
    for (int i=0; i<100; i++) {
        trial1 = randg.uniformShell(1.0,2.0);
        trial2 = randg.uniformShell(2.0,5.0);
        trial3 = randg.uniformShell(3.0,6.0);
        REQUIRE(trial1.norm() <= 2.0);
        REQUIRE(trial1.norm() >= 1.0);
        REQUIRE(trial2.norm() <= 5.0);
        REQUIRE(trial2.norm() >= 2.0);
        REQUIRE(trial3.norm() <= 6.0);
        REQUIRE(trial3.norm() >= 3.0);
    }
}

TEST_CASE("Particle class basic functionality", "[particle]") {
    double D = 1.0;
    double Drot = 0.5;
    auto position = vec3<double> {0.0, 0.0, 0.0};
    auto orientation = quaternion<double> {1.0, 0.0, 0.0, 0.0};
    particle part = particle(D, Drot, position, orientation);
    // Constructor consistency test
    REQUIRE(part.getD() == D);
    REQUIRE(part.getDrot() == Drot);
    REQUIRE(part.position == position);
    REQUIRE(part.orientation == orientation);
    // Some functionality testing
    auto newposition = vec3<double> {1.0, 2.5, 3.0};
    auto neworientation = quaternion<double> {0.0, 0.0, 0.0, 1.0};
    part.setNextPosition(newposition);
    part.setNextOrientation(neworientation);
    // Check main values haven't been changed
    REQUIRE(part.position == position);
    REQUIRE(part.orientation == orientation);
    // Update current values with new ones ("nextValues")
    part.updatePosition();
    part.updateOrientation();
    // Check values indeed were changed
    REQUIRE(part.position == newposition);
    REQUIRE(part.orientation == neworientation);
}

TEST_CASE("ParticleMS class basic functionality", "[particleMS]") {
    int type = 0;
    int state = 2;
    double D = 1.0;
    double Drot = 0.5;
    auto position = vec3<double> {0.0, 0.0, 0.0};
    auto orientation = quaternion<double> {1.0, 0.0, 0.0, 0.0};
    particleMS partMS = particleMS(type, state, D, Drot, position, orientation);
    // Constructor consistency test
    REQUIRE(partMS.getType() == type);
    REQUIRE(partMS.getState() == state);
    REQUIRE(partMS.getD() == D);
    REQUIRE(partMS.getDrot() == Drot);
    REQUIRE(partMS.position == position);
    REQUIRE(partMS.orientation == orientation);
    // Some functionality testing
    int newstate = 0;
    int newtype = 1;
    auto newposition = vec3<double> {1.0, 2.5, 3.0};
    auto neworientation = quaternion<double> {0.0, 0.0, 0.0, 1.0};
    partMS.setNextType(newtype);
    partMS.setNextState(newstate);
    partMS.setNextPosition(newposition);
    partMS.setNextOrientation(neworientation);
    // Check main values haven't been changed
    REQUIRE(partMS.getType() == type);
    REQUIRE(partMS.getState() == state);
    REQUIRE(partMS.position == position);
    REQUIRE(partMS.orientation == orientation);
    // Update current values with new ones ("nextValues")
    partMS.updatePosition();
    partMS.updateOrientation();
    partMS.updateType();
    partMS.updateState();
    // Check values indeed were changed
    REQUIRE(partMS.getType() == newtype);
    REQUIRE(partMS.getState() == newstate);
    REQUIRE(partMS.position == newposition);
    REQUIRE(partMS.orientation == neworientation);
}


TEST_CASE("Fundamental CTMSM parameters and propagation test", "[ctmsm]") {
    int msmid = 0;
    std::vector<std::vector<double>> tmatrix ={ {-5, 2, 3,}, {3, -6, 3}, {1, 3, -4} };
    // sum of outgoing rates for each state (row sum without diagonal negative value)
    std::vector<double> lambda0 = {5, 6, 4};
    // cumulative sum of outgoing rates for each state
    std::vector<std::vector<double>> ratescumsum = {{2, 5}, {3, 6}, {1, 4}};
    long seed = 0;
    ctmsm ctmsmTest = ctmsm(msmid, tmatrix, seed);
    REQUIRE(ctmsmTest.nstates == 3);
    REQUIRE(ctmsmTest.getLambda0() == lambda0);
    REQUIRE(ctmsmTest.getRatescumsum() == ratescumsum);
    // Propagation test requires creating a particleMS
    int type = 0;
    int state = 2;
    double D = 1.0;
    double Drot = 0.5;
    auto position = vec3<double> {0.0, 0.0, 0.0};
    auto orientation = quaternion<double> {1.0, 0.0, 0.0, 0.0};
    particleMS partMS = particleMS(type, state, D, Drot, position, orientation);
    // Create a second particle and a second ctmsm with same seed
    ctmsm ctmsmTest2 = ctmsm(msmid, tmatrix, seed);
    particleMS partMS2 = particleMS(type, state, D, Drot, position, orientation);
    /* Propagate each particle using the same seed but the update and noupdate method,
     * respectively, and compare output at each timestep. */
    for (int i=0; i<100; i++) {
        ctmsmTest.propagate(partMS, 3);
        ctmsmTest2.propagateNoUpdate(partMS2, 3);
        partMS2.updateState();
        REQUIRE(partMS.getState() == partMS2.getState());
    }
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

TEST_CASE("Write to trajectory", "[writeTraj]"){
    auto p1 = vec3<double> {0.0, 0.2, 0.0};
    auto p2 = vec3<double> {0.0, 0.4, 0.0};
    auto o1 = quaternion<double> {1.0, 0.0, 0.0, 0.0};
    auto o2 = quaternion<double> {1.0, 0.0, 0.0, 0.0};
    particle part1 = particle(1., 1., p1, o1);
    particle part2 = particle(1., 1., p2, o2);
    std::vector<particle> particles{part1, part2};
    trajectoryPositionOrientation traj(2, 1);
    traj.sample(0.1, particles); //traj.sample(0.1, particles);
}

TEST_CASE("Fundamental trajectory recording", "[trajectory]") {
    auto p1 = vec3<double> {0.0, 0.0, 0.0};
    auto p2 = vec3<double> {0.0, 0.0, 0.0};
    auto o1 = quaternion<double> {1.0, 0.0, 0.0, 0.0};
    auto o2 = quaternion<double> {1.0, 0.0, 0.0, 0.0};
    particle part1(1., 1., p1, o1);
    particle part2(1., 1., p2, o2);
    std::vector<particle> particles {part1, part2};
    overdampedLangevin integrator(0.01, 15, "rigidbody");
    simulation sim(integrator);
    std::string trajtype = "positionOrientation";
    sim.run(particles, 10000, 10, 1024, "test.h5", false, false, false, trajtype);
    auto data = sim.traj->getTrajectoryData();
    REQUIRE(data.size() == 2000);
}

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