//
// Created by maojrs on 8/16/18.
//

//#include <cmath>
#include "particle.hpp"
#include "vec3.hpp"
#include "potentials/gaussians3D.hpp"

namespace msmrd {
    /**
     * Class implementation for Gaussians3D external potential
     * @param nminima number of minimas/gaussians in the potential
     * @param maxrad maximum radius where all centers of Gaussians will be contained
     * @param scalefactor adjusts strength of the potential muliplying by a constant
     * @param seed variable for random number generation (Note seed = -1 corresponds to random device)
     * Constructor, sets random location and std. deviation of nminima Gaussians
     */
    gaussians3D::gaussians3D(unsigned long nminima, double maxrad, double scalefactor, long seed)
            : nminima(nminima), maxrad(maxrad), scalefactor(scalefactor), seed(seed) {
        randg.setSeed(seed);
        minimas.resize(nminima);
        stddevs.resize(nminima);
        for (int i = 0; i < nminima; i++) {
            minimas[i] = randg.uniformSphere(maxrad);
            stddevs[i][0] = randg.uniformRange(0.1 * maxrad, 0.3 * maxrad);
            stddevs[i][1] = randg.uniformRange(0.1 * maxrad, 0.3 * maxrad);
            stddevs[i][2] = randg.uniformRange(0.1 * maxrad, 0.3 * maxrad);
        }
    }


    // Returns value of potential at position x
    double gaussians3D::evaluate(const particle &part) {
        vec3<double> x = part.position;
        double output = 0;
        double gauss;
        double denom;
        vec3<double> m, sig;
        for (int i = 0; i < nminima; i++) {
            m = 1.0 * minimas[i];
            sig = 1.0 * stddevs[i];
            gauss = std::exp(-std::pow(x[0] - m[0], 2) / (2 * std::pow(sig[0], 2))
                             - std::pow(x[1] - m[1], 2) / (2 * std::pow(sig[1], 2))
                             - std::pow(x[2] - m[2], 2) / (2 * std::pow(sig[2], 2)));
            denom = std::pow(2 * M_PI, 3.0 / 2.0) * sig[0] * sig[1] * sig[2];
            output -= gauss / denom;
        }
        return scalefactor * output;
    };


    // Returns minus gradient of potential (force) at position x and zero torque
    std::array<vec3<double>, 2> gaussians3D::forceTorque(const particle &part) {
        vec3<double> x = part.position;
        vec3<double> force = vec3<double>(0, 0, 0);
        vec3<double> torque = vec3<double>(0, 0, 0);
        vec3<double> m, sig, grad;
        double gauss, denom;
        for (int i = 0; i < nminima; i++) {
            m = 1.0 * minimas[i];
            sig = 1.0 * stddevs[i];
            gauss = std::exp(-std::pow(x[0] - m[0], 2) / (2 * std::pow(sig[0], 2))
                             - std::pow(x[1] - m[1], 2) / (2 * std::pow(sig[1], 2))
                             - std::pow(x[2] - m[2], 2) / (2 * std::pow(sig[2], 2)));
            denom = std::pow(2 * M_PI, 3.0 / 2.0) * sig[0] * sig[1] * sig[2];
            grad[0] = -(2 * (x[0] - m[0]) / (2 * std::pow(sig[0], 2))) * gauss / denom;
            grad[1] = -(2 * (x[1] - m[1]) / (2 * std::pow(sig[1], 2))) * gauss / denom;
            grad[2] = -(2 * (x[2] - m[2]) / (2 * std::pow(sig[2], 2))) * gauss / denom;
            force -= grad;
        }
        return {force, torque};
    }
}