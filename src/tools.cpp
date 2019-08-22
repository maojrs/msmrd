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
        auto resultquat = quaternion<double>(p);
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

    // Calculates relative distance (p2-p1) of two vectors (p1, p2) in a periodic box, returns relative distance.
    vec3<double> distancePeriodicBox(vec3<double> p1, vec3<double> p2, vec3<double> edgeslength) {
        vec3<double> p1Periodic = 1.0*p1;
        // Loop over three coordinates (x,y,z)
        for (int i = 0; i < 3; i++){
            if ( (p2[i] - p1[i]) > 0.5*edgeslength[i]) {
                p1Periodic[i] += edgeslength[i];
            }
            if ( (p2[i] - p1[i]) < -0.5*edgeslength[i]) {
                p1Periodic[i] -= edgeslength[i];
            }
        }
        return p2 - p1Periodic;
    }

    /* Calculates relative distance (p2-p1) of two vectors (p1, p2) in a periodic box, returns tuple with
     * "virtual" position of p1 as first entry and relative distance as second entry */
    std::array<vec3<double>, 2> distancePeriodicBoxComplete(vec3<double> p1, vec3<double> p2,
                                                                       vec3<double> edgeslength) {
        vec3<double> p1Periodic = 1.0*p1;
        // Loop over three coordinates (x,y,z)
        for (int i = 0; i < 3; i++){
            if (p2[i] - p1[i] > 0.5*edgeslength[i]) {
                p1Periodic[i] += edgeslength[i];
            }
            if (p2[i] - p1[i] < -0.5*edgeslength[i]) {
                p1Periodic[i] -= edgeslength[i];
            }
        }
        return {p1Periodic, p2 - p1Periodic};
    }

    // Calculates relative position between two particles taking into account possible periodic boundary.
    vec3<double> calculateRelativePosition(vec3<double> pos1, vec3<double> pos2, bool boundaryActive,
                                           std::string boundaryType, vec3<double> boxsize){
        vec3<double> relPosition;
        if (boundaryActive and boundaryType == "periodic") {
            relPosition = distancePeriodicBox(pos1, pos2, boxsize);
        } else {
            relPosition = pos2 - pos1;
        }
        return relPosition;
    }

    double dotQuaternion(quaternion<double> q0, quaternion<double> q1) {
        return q0[0]*q1[0] + q0[1]*q1[1] + q0[2]*q1[2]+ q0[3]*q1[3];
    }

    // Calculates slerp (spherical linear interpolation between two quaternions).
    quaternion<double> quaternionSlerp(quaternion<double> q0, quaternion<double> q1, double t) {
        const double dotThreshold = 0.9995;

        // Normalize since only unit quaternions are valid rotations.
        q0 = q0/q0.norm();
        q1 = q1/q1.norm();

        // Compute the cosine of the angle between the two vectors.
        double dot = dotQuaternion(q0, q1);

        /* If the dot product is negative, slerp won't take the shorter path. As v1 and -v1 are equivalent when
         * equivalent rotations, we only need reversing one quaternion. */
        if (dot < 0) {
            q1 = -1*q1;
            dot = -dot;
        }

        if (dot > dotThreshold) {
            // If the inputs are too close for comfort, linearly interpolate
            // and normalize the result.

            quaternion<double> result = q0 + t*(q1 - q0);
            result = result/result.norm();
            return result;
        }

        // Since dot is in range [0, dotThreshold], acos is safe
        double theta_0 = acos(dot);        // theta_0 = angle between input vectors
        double theta = theta_0*t;          // theta = angle between v0 and result
        double sin_theta = sin(theta);     // compute this value only once
        double sin_theta_0 = sin(theta_0); // compute this value only once

        double s0 = cos(theta) - dot * sin_theta / sin_theta_0;  // == sin(theta_0 - theta) / sin(theta_0)
        double s1 = sin_theta / sin_theta_0;

        return (s0 * q0) + (s1 * q1);
    }
}
