//
// Created by maojrs on 3/28/19.
//

#include <catch2/catch.hpp>
#include "trajectories/trajectory.hpp"
#include "trajectories/trajectoryPosition.hpp"
#include "trajectories/trajectoryPositionOrientation.hpp"
#include "integrators/overdampedLangevin.hpp"
#include "simulation.hpp"

using namespace msmrd;

TEST_CASE("Write to trajectory", "[writeTraj]"){
    auto p1 = vec3<double> {0.0, 0.2, 0.0};
    auto p2 = vec3<double> {0.0, 0.4, 0.0};
    auto o1 = quaternion<double> {1.0, 0.0, 0.0, 0.0};
    auto o2 = quaternion<double> {1.0, 0.0, 0.0, 0.0};
    particle part1 = particle(1., 1., p1, o1);
    particle part2 = particle(1., 1., p2, o2);
    std::vector<particle> particles{part1, part2};
    trajectoryPositionOrientation traj(2, 1);
    traj.sample(0.1, particles);
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
