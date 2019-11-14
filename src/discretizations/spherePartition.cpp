//
// Created by maojrs on 1/31/19.
//
#include "discretizations/spherePartition.hpp"

namespace msmrd{

    /*
     * Constructor creates a spherical equal area partition on the surface of a sphere.
     * @param numSections the number of sections that the sphere should be partioned in.
     */
    spherePartition::spherePartition(int numSections): numSections(numSections){
        partitionSphere();
    }

    // Calculates area of cap given polar angle.
    double spherePartition::angle2CapArea(double phi) {
        return (4 / scaling) * M_PI * (std::sin(phi / 2.0)) * (std::sin(phi / 2.0));
    }

    // Calculates polar angle corresponding to a given cap area
    double spherePartition::capArea2Angle(double area) {
        return 2.0 * std::asin(std::sqrt(area / (4.0 * M_PI / scaling)));
    }

    // Calculates state area (area of each region)
    double spherePartition::stateArea() {
        return (4*M_PI/numSections) / scaling;
    };

    /* Calculate equal area partition of unit sphere with "num_sections" sections in partition.
     *  Calculates three vectors. "regionsPerCollar" indicates number of sections in each
     *  collar (int). "phis" indicates location of cuts that define collars in the polar
     *  angle std::vector<double>. "thetas" denotes the location of azimutal cuts as for each collar,
     *  as a vector for each collar std::vector<std::vector<double>>>. (called by constructor)
     *  */
     void spherePartition::partitionSphere() {
        //Calculate areas of each state and polar caps angle (phi0 and pi-phi0)
        double state_area = stateArea();
        double phi0 = capArea2Angle(state_area);
        // Calculate the number of collars between the polar caps
        double ideal_collar_angle = std::sqrt(state_area);
        double ideal_num_collars = (M_PI - 2 * phi0) / ideal_collar_angle;
        int num_collars = static_cast<int>(std::max(1.0, std::round(ideal_num_collars)));
        if (numSections == 2) {
            num_collars = 0;
        }
        double collar_angle = (M_PI - 2 * phi0) / num_collars;
        // Initialize variables for number of regions in each collar
        std::vector<double> ideal_regionsPerCollar;
        ideal_regionsPerCollar.resize(num_collars);
        phis.resize(num_collars + 2);
        phis[0] = 0;
        regionsPerCollar.resize(num_collars);
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
            cap_area_phi1 = angle2CapArea(phi0 + i * collar_angle);
            cap_area_phi2 = angle2CapArea(phi0 + (i + 1) * collar_angle);
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
            phis[i + 1] = capArea2Angle(summ * state_area);
            phis[num_collars + 1] = M_PI - phi0;
            // Obtain list of thetas for a given collar
            regsPerCollar_i = static_cast<unsigned long>(regionsPerCollar[i]);
            std::vector<double> thetasi;
            thetasi.resize(regsPerCollar_i);
            dth = 2.0 * M_PI / regionsPerCollar[i] / scaling;
            for (int j = 0; j < regsPerCollar_i; j++) {
                thetasi[j] = j * dth;
            }
            thetas.push_back(thetasi);
        }
        regionsPerCollar.push_back(1);
        it = regionsPerCollar.begin();
        regionsPerCollar.insert(it, 1);
    }


    // Adds an offset on the thetas, must be a positive value and preferably smaller than all the possible dthetas
    void spherePartition::setThetasOffset(double offset) {
         double dtheta = std::abs(thetas[0][1] - thetas[0][0]);
         if ( offset < 0) {
             throw std::invalid_argument("Offset must be positive and preferably smaller than smaller "
                                         "dtheta in discretization");
         }
        thetasOffset = offset;
        for (auto &thetaList : thetas) {
             for (auto &theta : thetaList) {
                 theta += offset;
             }
         }
     };


    /* Assuming partitioned sphere sits in origin, given a vector coordinate, find section number
     * that corresponds to section that the line generated by the vector intersects. */
    int spherePartition::getSectionNumber(vec3<double> coordinate) {
        if (scaling == 2 and coordinate[1] < 0) {
            throw std::invalid_argument("Error: y coordinate must be positive in half sphere discretization");
        }
        // Calculate theta and phi of coordinate
        double theta = std::atan2(coordinate[1], coordinate[0]) - thetasOffset;
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
        if (currentCollarIndex == regionsPerCollar.size()-1) {
            sectionNum = numSections;
            return sectionNum;
        }
        std::vector<double> collarThetas = thetas[currentCollarIndex - 1];
        int numThetaCuts = collarThetas.size();
        for (int i = 0; i < numThetaCuts; i++){
            if(theta >= (collarThetas[numThetaCuts - 1 - i] - thetasOffset)){
                currentThetaIndex = numThetaCuts - 1 - i;
                break;
            }
        }
        sectionNum = std::accumulate(std::begin(regionsPerCollar), std::next(std::begin(regionsPerCollar),
                                     currentCollarIndex), 0) + currentThetaIndex + 1;

        return sectionNum;
    }

    /* Returns phi-angles (polar) and theta-angles (azimuthal) that correspond to the sectionnumber
     * in the sphere partition. Note if thetasOffset != 0, then it can return one thetas interval with
     * value larger than 2pi. However this should not affect execution of dependencies.*/
    std::tuple<std::array<double, 2>, std::array<double, 2>> spherePartition::getAngles(int secNumber) {
        if (secNumber > numSections) {
            throw std::invalid_argument("Error: section number is larger than number of partitions");
        }
        // Get collar
        int collar = 0;
        int sections = 1;
        while (sections < secNumber) {
            collar += 1;
            sections = std::accumulate(std::begin(regionsPerCollar),
                                       std::next(std::begin(regionsPerCollar), collar+1), 0);
        }
        // Find phis
        double phi1 = phis[collar];
        double phi2;
        if (collar + 1 < regionsPerCollar.size()) {
            phi2 = phis[collar + 1];
        } else {
            phi2 = M_PI;
        }
        // Find thetas
        double theta1;
        double theta2;
        std::vector<double> thetasCollar;
        int statesInCollar;
        int prevStates = std::accumulate(std::begin(regionsPerCollar),
                                         std::next(std::begin(regionsPerCollar), collar), 0);
        if ((prevStates == 0) or (prevStates == numSections - 1)) {
            theta1 = 0;
            theta2 = 2 * M_PI / scaling;
        } else {
            thetasCollar = thetas[collar - 1];
            statesInCollar = secNumber - prevStates;
            theta1 = thetasCollar[statesInCollar - 1];
            if (statesInCollar == thetasCollar.size()) {
                theta2 = (2 * M_PI + thetasOffset)/scaling;
            } else {
                theta2 = thetasCollar[statesInCollar];
            }
        }
        std::array<double, 2> phiInterval = {phi1, phi2};
        std::array<double, 2> thetaInterval = {theta1, theta2};
        return std::make_tuple(phiInterval, thetaInterval);
    }

    /*  Returns list of three vectors that define the partition and calculate with partionSphere functions.
     * The first vector indicates number of sections in each collar (int). The second one indicates location
     * of cuts that define collars in the polar angle std::vector<double>. The third one denotes the location
     * of azimutal cuts as for each collar, as a vector for each collar std::vector<std::vector<double>>> */
    std::tuple<std::vector<int>, std::vector<double>,
            std::vector<std::vector<double>>> spherePartition::getPartition() {
        return std::make_tuple(regionsPerCollar, phis, thetas);
    };

    // Pybind version of getSectionNumber, uses arbitrary vectors instead of vec3.
    int spherePartition::getSectionNumberPyBind(std::vector<double> coord) {
        return getSectionNumber(vec3<double>{coord[0], coord[1], coord[2]});
    }
}
