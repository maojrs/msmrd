//
// Created by maojrs on 10/2/18.
//

#pragma once
#include <vector>
#include "quaternion.hpp"
#include "vec3.hpp"

namespace msmrdtools {

    /*
     * Useful vector conversions.
     */

    /* From an array of N vectors (vec3) returns a vector of vectors (std::vector) with the same dimension. */
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


    /*
     * Useful quaternion functions.
     */


    template<typename scalar=double>
    quaternion<scalar> axisangle2quaternion(const vec3<double> &phi) {
        double phinorm = phi.norm();
        if (phinorm != 0) {
            vec3<double> phiunit = phi / phinorm;
            double s = cos(0.5 * phinorm);
            double p = sin(0.5 * phinorm);
            return {s, p * phiunit[0], p * phiunit[1], p * phiunit[2]};
        } else {
            //returns unit quaternion (no rotation)
            return {1, 0, 0, 0};
        }
    };

    template<typename scalar=double>
    vec3<double> rotateVec(vec3<double> p, quaternion<double> q) {
        vec3<double> result;
        quaternion<double> resultquat = 1.0*quaternion<double>(p);
        resultquat = q*(resultquat*q.conj());
        result[0] = resultquat[1];
        result[1] = resultquat[2];
        result[2] = resultquat[3];
        return result;
    };
}

