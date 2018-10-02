//
// Created by dibakma on 26.06.18.
//
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include <catch2/catch.hpp>
#include <functional>
#include "msm.hpp"
#include "integrators/odLangevin.hpp"
#include "simulation.hpp"
#include "trajectory.hpp"
#include "particle.hpp"
#include "quaternion.hpp"
#include "randomgen.hpp"
#include "vec3.hpp"

using namespace msmrd;
using msm = msmrd::discreteTimeMarkovStateModel;
using ctmsm = msmrd::continuousTimeMarkovStateModel;
using namespace std::placeholders;

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
    std::string bodytype = "rigidsolid";
    auto position = vec3<double> {0.0, 0.0, 0.0};
    auto orientation = quaternion<double> {1.0, 0.0, 0.0, 0.0};
    particle part = particle(D, Drot, bodytype, position, orientation);
    // Constructor consistency test
    REQUIRE(part.getD() == D);
    REQUIRE(part.getDrot() == Drot);
    REQUIRE(part.getBodyType() == bodytype);
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
    std::string bodytype = "rigidsolid";
    auto position = vec3<double> {0.0, 0.0, 0.0};
    auto orientation = quaternion<double> {1.0, 0.0, 0.0, 0.0};
    particleMS partMS = particleMS(type, state, D, Drot, bodytype, position, orientation);
    // Constructor consistency test
    REQUIRE(partMS.getType() == type);
    REQUIRE(partMS.getState() == state);
    REQUIRE(partMS.getD() == D);
    REQUIRE(partMS.getDrot() == Drot);
    REQUIRE(partMS.getBodyType() == bodytype);
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
    std::string bodytype = "rigidsolid";
    auto position = vec3<double> {0.0, 0.0, 0.0};
    auto orientation = quaternion<double> {1.0, 0.0, 0.0, 0.0};
    particleMS partMS = particleMS(type, state, D, Drot, bodytype, position, orientation);
    // Create a second particle and a second ctmsm with same seed
    ctmsm ctmsmTest2 = ctmsm(msmid, tmatrix, seed);
    particleMS partMS2 = particleMS(type, state, D, Drot, bodytype, position, orientation);
    /* Propagate each particle using the same seed but the update and noupdate method,
     * respectively, and compare output at each timestep. */
    for (int i=0; i<100; i++) {
        ctmsmTest.propagate(partMS, 3);
        ctmsmTest2.propagateNoUpdate(partMS2, 3);
        partMS2.updateState();
        REQUIRE(partMS.getState() == partMS2.getState());
    }
}

TEST_CASE("Write to trajectory", "[writeTraj]"){
    auto p1 = vec3<double> {0.0, 0.2, 0.0};
    auto p2 = vec3<double> {0.0, 0.4, 0.0};
    auto o1 = quaternion<double> {1.0, 0.0, 0.0, 0.0};
    auto o2 = quaternion<double> {1.0, 0.0, 0.0, 0.0};
    particle part1 = particle(1., 1., "rigidsolid", p1, o1);
    particle part2 = particle(1., 1., "rigidsolid", p2, o2);
    std::vector<particle> particles{part1, part2};
    trajectoryPositionOrientation traj(2, 1);
    traj.sample(0.1, particles); //traj.sample(0.1, particles);
}

TEST_CASE("Fundamental trajectory recording", "[trajectory]") {
    auto p1 = vec3<double> {0.0, 0.0, 0.0};
    auto p2 = vec3<double> {0.0, 0.0, 0.0};
    auto o1 = quaternion<double> {1.0, 0.0, 0.0, 0.0};
    auto o2 = quaternion<double> {1.0, 0.0, 0.0, 0.0};
    particle part1(1., 1., "rigidsolid", p1, o1);
    particle part2(1., 1., "rigidsolid", p2, o2);
    std::vector<particle> particles {part1, part2};
    odLangevin integrator(0.01, 15, true);
    simulation sim(integrator, particles);
    trajectoryPositionOrientation traj(2, 20000);
    sim.run(10000, traj, 1);
    REQUIRE(traj.data.size() == 20000);
}