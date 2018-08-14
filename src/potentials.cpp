//
// Created by maojrs on 8/14/18.
//
#include "randomgen.hpp"
#include "vec3.hpp"
#include "potentials.hpp"


// Constructor, sets random location and width of nminima Gaussians
gaussians3D::gaussians3D(int nminima, double maxrad, double scalfactor, long seed)
        : nminima(nminima), maxrad(maxrad), scalefactor(scalefactor), seed(seed) {
    randg.setSeed(seed);
    minimas.resize(nminima);
    stddevs.resize(nminima);
    for (int i=0; i<nminima; i++) {
        minimas[i] = randg.uniformSphere(maxrad);
        stddevs[i][0] = randg.uniformRange(0.1*maxrad,0.3*maxrad);
        stddevs[i][1] = randg.uniformRange(0.1*maxrad,0.3*maxrad);
        stddevs[i][2] = randg.uniformRange(0.1*maxrad,0.3*maxrad);
    }
}

// Returns value of potential at position x
double gaussians3D::evaluate(vec3<double> x) {
    double gauss = 0;
    double output = 0;
    for (int i=0; i<nminima; i++) {
        gauss = std::exp( -std::pow(x[0] - minimas[i][0],2)/(2*std::pow(stddevs[i][0],2)));
        gauss += std::exp(-std::pow(x[1] - minimas[i][1],2)/(2*std::pow(stddevs[i][1],2)));
        gauss += std::exp(-std::pow(x[2] - minimas[i][2],2)/(2*std::pow(stddevs[i][2],2)));
        output -= gauss;
    }
    return scalefactor*output;
};

// Needed for PyBinding evaluate, since it can take vectors as input.
double gaussians3D::evaluatePyBind(std::vector<double> pos) {
    vec3<double> x = vec3<double>(pos);
    return evaluate(x);
}


// Returns -gradient of potential (force) at position x
vec3<double> gaussians3D::force(vec3<double> x) {
    vec3<double> force = vec3<double>(0, 0, 0);
    vec3<double> m, sig, grad;
    double expall, denom;
    for (int i=0; i<nminima; i++) {
        m = 1.0*minimas[i];
        sig = 1.0*stddevs[i];
        expall = std::exp(-std::pow(x[0]-m[0],2)/(2*std::pow(sig[0],2))
                          -std::pow(x[1]-m[1],2)/(2*std::pow(sig[1],2))
                          -std::pow(x[2]-m[2],2)/(2*std::pow(sig[2],2)));
        denom = std::pow(2*M_PI, 3.0/2.0)*sig[0]*sig[1]*sig[2];
        grad[0] = -(2*(x[0]-m[0])/(2*std::pow(sig[0],2)))*expall/denom;
        grad[1] = -(2*(x[1]-m[1])/(2*std::pow(sig[1],2)))*expall/denom;
        grad[2] = -(2*(x[2]-m[2])/(2*std::pow(sig[2],2)))*expall/denom;
        force -= grad;
    }
    return scalefactor*force;
};

// Needed for PyBinding force, since it can take vectors as input.
std::vector<double> gaussians3D::forcePyBind(std::vector<double> pos) {
    vec3<double> x = vec3<double>(pos);
    vec3<double> forcex = force(x);
    std::vector<double> output;
    output[0] = forcex[0];
    output[1] = forcex[1];
    output[2] = forcex[2];
    return output;
}
