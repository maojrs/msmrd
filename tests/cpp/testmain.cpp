//
// Created by dibakma on 26.06.18.
//
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include <catch2/catch.hpp>
#include "simulation.hpp"
#include "vec3.hpp"
#include "particle.hpp"

unsigned int Factorial( unsigned int number ) {
    return number <= 1 ? number : Factorial(number-1)*number;
}

TEST_CASE( "Factorials are computed", "[factorial]" ) {
REQUIRE( Factorial(1) == 1 );
REQUIRE( Factorial(2) == 2 );
REQUIRE( Factorial(3) == 6 );
REQUIRE( Factorial(10) == 3628800 );
}

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