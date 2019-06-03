//
// Created by maojrs on 5/29/19.
//

#include <catch2/catch.hpp>
#include "particle.hpp"
#include "boundaries/boundary.hpp"
#include "boundaries/noBoundary.hpp"
#include "vec3.hpp"

using namespace msmrd;


TEST_CASE("No boundary class test", "[noBoundary]") {
    boundary *domainBoundary;
    auto defaultBoundary = noBoundary();
    domainBoundary = &defaultBoundary;
    REQUIRE( domainBoundary->getBoundaryType() == "none" );
    REQUIRE( domainBoundary->boxsize == vec3<double>(0., 0., 0.) );
}
