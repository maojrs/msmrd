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


    /*
    * C++ version of fucntions to create sphere partition of equal area. It is a c++ copy of the python
    * code in module msmrd2.tools.spherePartition.
    */
    namespace spherePartition {

        // Calculates area of cap given polar angle
        double angle_to_cap_area(double phi) {
            return 4 * M_PI * (std::sin(phi / 2.0)) * (std::sin(phi / 2.0));
        }

        // Calculates polar angle corresponding to a given cap area
        double cap_area_to_angle(double area) {
            return 2.0 * std::asin(std::sqrt(area / (4.0 * M_PI)));
        }

        /* Calculate equal area partition of unit sphere with "num_partitions" partitions.
         *  Returns list of three vectors. The first vector indicates number of sections in each
         *  collar (int). The second one indicates location of cuts that define collars in the polar
         *  angle std::vector<double>. The third one denotes the location of azimutal cuts as for each collar,
         *  as a vector for each collar std::vector<std::vector<double>>> .
         *  */
        std::tuple<std::vector<int>, std::vector<double>,
                std::vector<std::vector<double>>> partitionSphere(int num_partitions) {
            //Calculate areas of each state and polar caps angle (phi0 and pi-phi0)
            double state_area = 4 * M_PI / num_partitions;
            double phi0 = cap_area_to_angle(state_area);
            // Calculate the number of collars between the polar caps
            double ideal_collar_angle = std::sqrt(state_area);
            double ideal_num_collars = (M_PI - 2 * phi0) / ideal_collar_angle;
            int num_collars = static_cast<int>(std::max(1.0, std::round(ideal_num_collars)));
            if (num_partitions == 2) {
                num_collars = 0;
            }
            double collar_angle = (M_PI - 2 * phi0) / num_collars;
            // Initialize variables for number of regions in each collar
            std::vector<double> ideal_regionsPerCollar;
            ideal_regionsPerCollar.resize(num_collars);
            std::vector<double> phis;
            phis.resize(num_collars + 2);
            phis[0] = 0;
            std::vector<int> regionsPerCollar;
            regionsPerCollar.resize(num_collars);
            std::vector<std::vector<double>> thetas;
            thetas.resize(0);
            std::vector<double> a{0};
            /* Iterate over each collar to get right number of regions per collar
             * and correct location of phi angles of each collar. */
            double cap_area_phi1;
            double cap_area_phi2;
            double suma;
            double summ;
            unsigned long regsPerCollar_i;
            double dth;
            std::vector<int>::iterator it;
            for (int i = 0; i < num_collars; i++) {
                // Calculate num of regions in collar i
                cap_area_phi1 = angle_to_cap_area(phi0 + i * collar_angle);
                cap_area_phi2 = angle_to_cap_area(phi0 + (i + 1) * collar_angle);
                ideal_regionsPerCollar[i] = (cap_area_phi2 - cap_area_phi1) / state_area;
                regionsPerCollar[i] = static_cast<int>(std::round(ideal_regionsPerCollar[i] + a[i]));
                // Correct values of phi around collar i
                suma = 0;
                for (int j = 0; j < i + 1; j++) {
                    suma = suma + ideal_regionsPerCollar[j] - regionsPerCollar[j];
                }
                a.push_back(suma);
                summ = 1;
                for (int j = 0; j < i; j++) {
                    summ = summ + regionsPerCollar[j];
                }
                phis[i + 1] = cap_area_to_angle(summ * state_area);
                phis[num_collars + 1] = M_PI - phi0;
                // Obtain list of thetas for a given collar
                regsPerCollar_i = static_cast<unsigned long>(regionsPerCollar[i]);
                std::vector<double> thetasi;
                thetasi.resize(regsPerCollar_i);
                dth = 2.0 * M_PI / regionsPerCollar[i];
                for (int j = 0; j < regsPerCollar_i; j++) {
                    thetasi[j] = j * dth;
                }
                thetas.push_back(thetasi);
            }
            regionsPerCollar.push_back(1);
            it = regionsPerCollar.begin();
            regionsPerCollar.insert(it, 1);
            // return number of regions for all collars,
            // phi angles of collars and theta angles for each collar
            return std::make_tuple(regionsPerCollar, phis, thetas);
        }

        /* Assuming partioned sphere sits in origin, given a vector coordinate, find section number
         * that corresponds to section that the line generated by the vector intersects. */
        int getSectionNumber(vec3<double> coordinate, int numPartitions ) {
            // Get sphere partition
            auto spherePartition = msmrdtools::spherePartition::partitionSphere(numPartitions);
            auto numRegionsCollar = std::get<0>(spherePartition);
            auto phis = std::get<1>(spherePartition);
            auto thetas = std::get<2>(spherePartition);
            // Calculate theta and phi of coordinate
            double theta = std::atan2(coordinate[1], coordinate[0]);
            if (theta < 0) {
                theta += 2 * M_PI;
            }
            double r = coordinate.norm();
            double phi = std::acos(coordinate[2] / r);
            int currentCollarIndex;
            int currentThetaIndex;
            // Find intersection of coordinate with section
            int sectionNum;
            int numCollars = phis.size();
            for (int i = 0; i < numCollars; i++){
                if (phi >= phis[numCollars - 1 - i]){
                    currentCollarIndex = numCollars-1-i;
                    break;
                }
            }
            if (currentCollarIndex == 0) {
                sectionNum = 1;
                return sectionNum;
            }
            if (currentCollarIndex == numRegionsCollar.size()-1) {
                sectionNum = numPartitions;
                return sectionNum;
            }
            std::vector<double> collarThetas = thetas[currentCollarIndex - 1];
            int numThetaCuts = collarThetas.size();
            for (int i = 0; i < numThetaCuts; i++){
                if(theta >= collarThetas[numThetaCuts - 1 - i]){
                    currentThetaIndex = numThetaCuts - i - 1;
                    break;
                }
            }
            sectionNum = std::accumulate(std::begin(numRegionsCollar),
                    std::next(std::begin(numRegionsCollar), currentCollarIndex), 0) + currentThetaIndex + 1;

            return sectionNum;
        }

//                def
//        getAngles(secNumber, numPartitions = None
//        ):
//        '''
//        Return angles[phi1, phi2]
//        and [theta1, theta2]
//        of a
//        given section
//        number
//                in
//        the spherical
//        partition
//        '''
//        if numPartitions == None:
//        numPartitions == 20
//        if secNumber > numPartitions:
//        print("Error: section number is larger than number of partitions")
//        return
//        numRegionsCollar, phis,
//        thetas = partitionSphere(numPartitions)
//        collar = 0
//        while (
//        sum(numRegionsCollar[0
//        :collar+1])) < secNumber:
//        collar += 1
//        phi1 = phis[collar]
//        if collar+1 <
//        len(numRegionsCollar)
//        :
//        phi2 = phis[collar + 1]
//        else:
//        phi2 = np.pi
//# Find thetas
//        prevStates = sum(numRegionsCollar
//        [0:collar])
//        if ((prevStates == 0) or (prevStates == numPartitions -1)):
//        theta1 = 0
//        theta2 = 2 * np.pi
//        else:
//        thetasCollar = thetas[collar - 1]
//        statesInCollar = secNumber - prevStates
//        theta1 = thetasCollar[statesInCollar - 1]
//        if statesInCollar ==
//        len(thetasCollar)
//        :
//        theta2 = 2 * np.pi
//        else:
//        theta2 = thetasCollar[statesInCollar]
//        phiInterval = [phi1, phi2]
//        thetaInterval = [theta1, theta2]
//        return phiInterval, thetaInterval
//
    }

}
