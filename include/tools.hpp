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
    template <unsigned long N>
    double stdvecNorm(std::array<double, N> a, std::array<double, N> b) {
        if (a.size() != b.size()) {
            std::range_error(" Vectors need to be the same size to obtain a valid norm");
        }
        int vecsize = a.size();
        double diff = 0;
        for (int i = 0; i < vecsize; i++) {
            diff += (a[i] - b[i])*(a[i] - b[i]);
        }
        diff = std::sqrt(diff);
        return diff;
    };


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
    std::array<vec3<double>, 2> distancePeriodicBoxComplete(vec3<double> p1, vec3<double> p2,
                                                                       vec3<double> edgeslength);

    // Calculates relative position between two particles taking into account possible periodic boundary
    vec3<double> calculateRelativePosition(vec3<double> pos1, vec3<double> pos2, bool boundaryActive,
                                           std::string boundaryType, vec3<double> boxsize);

    // Calculates dot product between quaternions
    double dotQuaternion(quaternion<double> q0, quaternion<double> q1);

    // Calculates slerp (spherical linear interpolation) between quaternions, t between (0,1).
    quaternion<double> quaternionSlerp(quaternion<double> q0, quaternion<double> q1, double t);

}

