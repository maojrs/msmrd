//
// Created by maojrs on 10/2/18.
//

#pragma once
#include <vector>
#include <tuple>
#include "quaternion.hpp"
#include "vec3.hpp"

namespace msmrdtools {

    /*
     * Useful vector tools.
     */

    /* From an array of N vectors (vec3) returns a vector of vectors (std::vector) with the same dimension.
     * (template must be implemented in header.) */
    template<long unsigned int N>
    std::vector<std::vector<double>> array2Dtovec2D(std::array<vec3<double>, N> forceTorque) {
        std::vector<std::vector<double>> output;
        output.resize(N);
        for (int i = 0; i < N; i++) {
            output[i].resize(3);
            output[i][0] = 1.0*forceTorque[i][0];
            output[i][1] = 1.0*forceTorque[i][1];
            output[i][2] = 1.0*forceTorque[i][2];
        }
        return output;
    }

    // Calculates norm between two vectors of the same size
    double stdvecNorm(std::vector<double> a, std::vector<double> b);


    /*
     * Useful quaternion functions definitions.
     */

    /* Converts vector defining axis of rotation with length equal
     * to angle of rotation to its quaternion representation. */
    quaternion<double> axisangle2quaternion(const vec3<double> &phi);

    // Converts quaternion to axis-angle representation.
    vec3<double> quaternion2axisangle(const quaternion<double> q);

    // Rotates vector p by rotation represented by quaternion q.
    vec3<double> rotateVec(vec3<double> p, quaternion<double> q);

    // Calculate distance between quaternions/rotations
    double quaternionDistance(quaternion<double> q1, quaternion<double> q2);

    // Calculate minimum rotation angle along some axis to reach q2 from q1.
    double quaternionAngleDistance(quaternion<double> q1, quaternion<double> q2);

    // Calculates relative distance of two vectors (p1, p2) in a periodic box
    vec3<double> distancePeriodicBox(vec3<double> p1, vec3<double> p2, vec3<double> edgeslength);

    /* Calculates relative distance of two vectors (p1, p2) in a periodic box and returns
     * virtual position of p1 + relative distance */
    std::tuple<vec3<double>, vec3<double>> distancePeriodicBoxComplete(vec3<double> p1, vec3<double> p2,
                                                                       vec3<double> edgeslength);

}

