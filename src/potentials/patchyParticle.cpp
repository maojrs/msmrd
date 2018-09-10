//
// Created by maojrs on 9/10/18.
//

#include "potentials/patchyParticle.hpp"




namespace msmrd{ ;
    /*
     * @param radius radius at which patches are placed. Corresponds to the radius of the particle
     * @param patchesCoordintates list of patches coordinates in a sphere of unit radius
     */
    patchyParticle::patchyParticle(double sigma, std::vector<vec3<double>> patchesCoordinates)
            : sigma(sigma), patchesCoordinates(patchesCoordinates)  {}

    patchyParticle::patchyParticle(double sigma, std::vector<std::vector<double>> patchesCoords)
            : sigma(sigma) {
        for (int i=0; i < patchesCoords.size(); i++ ) {
            patchesCoordinates[i] = vec3<double> (patchesCoords[i]);
        }
    }


    double patchyParticle::evaluate(vec3<double> pos1, vec3<double> pos2, quaternion<double> u1, quaternion<double> u2) {}
    std::array<vec3<double>, 2> patchyParticle::forceTorque(vec3<double> pos1, vec3<double> pos2, quaternion<double> u1, quaternion<double> u2) {}

    /* Custom quadratic potential:
     * @param r is distance between particles,
     * @param eps strength of the potential
     * @param rstar distance to switch to second potential function
     * @param a the stiffness, together with rstar determines the range of the potential */
    double patchyParticle::quadraticPotential(double r, double eps, double rstar, double a) {
        // Parameters to force continuity and continuous derivative
        double rcritical = std::pow(sigma, 2)/(a*rstar);
        double b = (1.0 - a*std::pow(rstar/sigma, 2))/std::pow(rcritical/sigma - rstar/sigma, 2);
        // Calculate different regions of potential
        if (r < 0) {
            return eps;
        } else if (r < rstar) {
            return eps * (1.0 - a * std::pow((r / sigma),2));
        } else if (r < rcritical) {
            return eps * b * std::pow((rcritical / sigma - r / sigma), 2);
        } else {
            return 0.0;
        }
    }

    double patchyParticle::derivativeQuadraticPotential(double r, double rstar, double eps, double a) {
        // Parameters to force continuity and continuous derivative
        double rcritical = std::pow(sigma, 2)/(a*rstar);
        double b = (1.0 - a*std::pow(rstar/sigma, 2))/std::pow(rcritical/sigma - rstar/sigma, 2);
        // Calculate different regions of potential
        if (r < 0) {
            return 0.0;
        } else if (r < rstar) {
            return - 2.0*eps*a * r/std::pow(sigma,2);
        } else if (r < rcritical) {
            return -2.0*eps * b * (rcritical - r)/std::pow(sigma,2);
        } else {
            return 0.0;
        }
    }


}

