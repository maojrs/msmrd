//
// Created by maojrs on 1/29/19.
//

#include "tools.hpp"

namespace msmrdtools {

    /*
     * Implementation of tools functions defined in tools.hpp
     */

    // Calculates norm between two vectors of the same size
    double stdvecNorm(std::vector<double> a, std::vector<double> b) {
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
    }

    /* Converts vector defining axis of rotation with length equal
     * to angle of rotation to its quaternion representation. */
    quaternion<double> axisangle2quaternion(const vec3<double> &phi) {
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
    }

    /* Converts quaternion to angle-axis representation, where
     * the vector indicates axis of rotation and length indicate degrees of rotation. */
    vec3<double> quaternion2axisangle(const quaternion<double> q) {
        vec3<double> phi = {q[1], q[2], q[3]};
        double qnorm = phi.norm();
        phi = phi / qnorm;
        double theta = 2 * std::atan2(qnorm, q[0]);
        phi = phi * theta;
        return phi;
    }

    // Rotates vector p by rotation represented by quaternion q.
    vec3<double> rotateVec(vec3<double> p, quaternion<double> q) {
        vec3<double> result;
        quaternion<double> resultquat = 1.0*quaternion<double>(p);
        resultquat = q*(resultquat*q.conj());
        result[0] = resultquat[1];
        result[1] = resultquat[2];
        result[2] = resultquat[3];
        return result;
    }

    /* Calculate distance between quaternions/rotations. Gives 0 when rotations are the same
     * and 1 when represent rotations 180 degress apart. */
    double quaternionDistance(quaternion<double> q1, quaternion<double> q2) {
        double innerProduct = q1[0]*q2[0] + q1[1]*q2[1] + q1[2]*q2[2] + q1[3]*q2[3];
        return 1.0 - innerProduct*innerProduct;
    }

    /* Calculate minimum rotation angle along some axis to reach q2 from q1.
     * Output between 0 and 2pi. */
    double quaternionAngleDistance(quaternion<double> q1, quaternion<double> q2) {
        quaternion<double> relquat = q2 * q1.conj();
        vec3<double> relangle = quaternion2axisangle(relquat);
        return relangle.norm();
    }

    // Calculates relative distance (p2-p1) of two vectors (p1, p2) in a periodic box
    vec3<double> distancePeriodicBox(vec3<double> p1, vec3<double> p2, vec3<double> edgeslength) {
        vec3<double> p1Periodic = p1;
        vec3<double> halfedgelength = 0.5*edgeslength;
        // Loop over three coordinates (x,y,z)
        for (int i = 0; i < 3; i++){
            if (p2[i] - p1[i] > halfedgelength[i]) {
                p1Periodic[i] += edgeslength[i];
            }
            if (p2[i] - p1[i] < -halfedgelength[i]) {
                p1Periodic[i] -= edgeslength[i];
            }
        }
        return p2 - p1Periodic;
    }
}
