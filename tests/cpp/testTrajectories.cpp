//
// Created by maojrs on 3/28/19.
//

#include <catch2/catch.hpp>
#include "trajectories/trajectory.hpp"
#include "trajectories/trajectoryPosition.hpp"
#include "trajectories/trajectoryPositionOrientation.hpp"
#include "trajectories/discrete/patchyProteinTrajectory.hpp"
#include "integrators/overdampedLangevin.hpp"
#include "simulation.hpp"
#include "tools.hpp"


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

TEST_CASE("Patchy Protein trajectory", "[patchyProteinTrajectory]") {
    /* Define relative position vectors measured from particle 1, as in setBoundStates()*/
    std::array<vec3<double>, 6> relPos;
    relPos[0] = {1., 0., 0.};
    relPos[1] = {0., 1., 0.};
    relPos[2] = {0., 0., 1.};
    relPos[3] = {-1., 0., 0.};
    relPos[4] = {0., -1., 0.};
    relPos[5] = {0., 0., -1.};
    /* Relative rotations (assuming particle 1 fixed) of particle 2 that yield the 6 bound states
     * in the axis-angle representation. (One needs to make drawing to understand), as defined in setBoundStates()*/
    std::array<vec3<double>, 6> rotations;
    rotations[0] = {0.0, 0.0, M_PI}; //ok
    rotations[1] = {0.0, 0.0, -M_PI / 2.0}; //ok
    rotations[2] = {0.0, M_PI / 2.0, 0.0}; //ok
    rotations[3] = {0.0, 0.0, 0.0}; //ok
    rotations[4] = {0.0, 0.0, M_PI / 2.0}; //ok
    rotations[5] = {0.0, -M_PI / 2.0, 0.0}; //ok
    /*Convert rotations in the axis angle representation to quaternions */
    std::array<quaternion<double>, 6> quatRotations;
    for (int i = 0; i < 6; i++) {
        quatRotations[i] = msmrdtools::axisangle2quaternion(rotations[i]);
    }
    // Define patchy protein trajectory
    patchyProteinTrajectory traj(2,1);
    // Check states calculated by sampleDiscreteState fucntion match the states defined originally in setBoundStates()
    for (int i=0; i<6; i ++){
        auto p1 = vec3<double> {0.0, 0.0, 0.0};
        auto p2 = vec3<double> {relPos[i]};
        auto o1 = quaternion<double> {1.0, 0.0, 0.0, 0.0};
        auto o2 = quaternion<double> {quatRotations[i]};
        particle part1(1., 1., p1, o1);
        particle part2(1., 1., p2, o2);
        auto discreteState = traj.sampleDiscreteState(part1,part2);
        REQUIRE(discreteState == i+1);
    }
}
