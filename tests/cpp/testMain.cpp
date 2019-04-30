//
// Created by dibakma on 26.06.18.
//
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include <catch2/catch.hpp>
#include "eventManager.hpp"
#include "particle.hpp"
#include "quaternion.hpp"
#include "randomgen.hpp"
#include "simulation.hpp"
#include "tools.hpp"
#include "vec3.hpp"

using namespace msmrd;

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

TEST_CASE("Quaternion Slerp", "[quaternions]") {
    double t = 0.7; //between 0 and 1
    double tolerance = 0.00001;
    quaternion<double> q1({1.,2.,3.,4.});
    quaternion<double> q2({5.,1.,7.,5.});
    q1 = q1/q1.norm();
    q2 = q2/q2.norm();
    quaternion<double> qslerp1 = msmrdtools::quaternionSlerp(q1, q2, t);
    quaternion<double> qslerp2 = msmrdtools::quaternionSlerp(q2, q1, 1 - t);
    auto qerror = qslerp1 - qslerp2;
    REQUIRE( qerror.norm() < tolerance );
    auto qrel1 = qslerp1*q1.conj();
    auto qrel2 = qslerp2*q2.conj();
    auto qerror2 = qrel1*q1 - qrel2*q2;
    REQUIRE( qerror2.norm() < tolerance  );
//  REQUIRE( q1*q2 == v4 );
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

TEST_CASE("Event manager functionality", "[eventManager]") {
    eventManager eventMgr = eventManager();
    double waitTime = 5.5;
    int endState = 2;
    std::string inORout = "in";
    // Test adding events
    eventMgr.addEvent(waitTime, endState, 1, 2, inORout);
    eventMgr.addEvent(0.5*waitTime, 2*endState, 3, 5, inORout);
    eventMgr.addEvent(2*waitTime, endState/2, 7, 6, inORout);
    eventMgr.addEvent(3*waitTime, 3*endState, 3, 1, inORout);
    REQUIRE(eventMgr.eventList.size() == 4);
    // Test sorting of events by waiting time
    eventMgr.sort();
    REQUIRE(eventMgr.getEventTime(0) == 0.5*waitTime);
    REQUIRE(eventMgr.getEventTime(1) == waitTime);
    REQUIRE(eventMgr.getEventTime(2) == 2*waitTime);
    REQUIRE(eventMgr.getEventTime(3) == 3*waitTime);
    // Test getting and removing events
    auto someEvent = eventMgr.getEvent(2);
    eventMgr.removeEvent(someEvent);
    REQUIRE(eventMgr.eventList.size() == 3);
    REQUIRE(eventMgr.getEventTime(0) == 0.5*waitTime);
    REQUIRE(eventMgr.getEventTime(1) == waitTime);
    REQUIRE(eventMgr.getEventTime(2) == 3*waitTime);
    eventMgr.removeEvent(1);
    REQUIRE(eventMgr.eventList.size() == 2);
    REQUIRE(eventMgr.getEventTime(0) == 0.5*waitTime);
    REQUIRE(eventMgr.getEventTime(1) == 3*waitTime);
    // Test advance time
    eventMgr.addEvent(waitTime, endState, 1, 2, inORout);
    eventMgr.addEvent(2*waitTime, endState/2, 7, 6, inORout);
    eventMgr.sort();
    eventMgr.advanceTime(1.5);
    REQUIRE(eventMgr.getEventTime(0) == 0.5*waitTime - 1.5);
    REQUIRE(eventMgr.getEventTime(1) == waitTime - 1.5);
    REQUIRE(eventMgr.getEventTime(2) == 2*waitTime - 1.5);
    REQUIRE(eventMgr.getEventTime(3) == 3*waitTime - 1.5);

}
