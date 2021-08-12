//
// Created by maojrs on 8/12/21.
//


//#include <cmath>
#include "particle.hpp"
#include "vec3.hpp"
#include "potentials/harmonic.hpp"

namespace msmrd {
    /**
     * Class implementation for harmonic external potential
     */

    /* Alternative constructor to specify location of minimas and standard deviation explicitly. */
    harmonic::harmonic(std::vector<double> minima_in, std::vector<double> kconstants_in, double scalefactor)
    : scalefactor(scalefactor) {
        if (minima_in.size() !=  3 or kconstants_in.size() !=3 ){
            throw std::invalid_argument("Minima and kconstants should be three-dimensional arrays");
        }
        minima = vec3<double>{minima_in[0], minima_in[1], minima_in[2]};
        kconstants = vec3<double>{kconstants_in[0], kconstants_in[1], kconstants_in[2]};
    }

    /* Alternative constructor to specify location of minimas and standard deviation explicitly, as well as
     * the specific particle types on which the potential will act on. */
    harmonic::harmonic(std::vector<double> minima_in, std::vector<double> kconstants_in,
            std::vector<int> partTypes, double scalefactor) :
            harmonic(minima_in, kconstants_in, scalefactor)  {
        particleTypes = partTypes;
    }


    // Returns value of potential at position x
    double harmonic::evaluate(particle &part) {
        vec3<double> x = part.position;
        vec3<double> m, k;
        double harmonicPot;
        if (particleTypes.empty() or std::find(particleTypes.begin(),
                                               particleTypes.end(), part.type) != particleTypes.end()) {
            m = 1.0 * minima;
            k = 1.0 * kconstants;
            harmonicPot = k[0] * std::pow(x[0] - m[0], 2) +
                          k[1] * std::pow(x[1] - m[1], 2) +
                          k[2] * std::pow(x[2] - m[2], 2) ;
        }
        return scalefactor * harmonicPot;
    };


    // Returns minus gradient of potential (force) at position x and zero torque
    std::array<vec3<double>, 2> harmonic::forceTorque(particle &part) {
        vec3<double> x = part.position;
        vec3<double> force = vec3<double>(0, 0, 0);
        vec3<double> torque = vec3<double>(0, 0, 0);
        vec3<double> m, k, grad;
        double harmonicPot;
        if (particleTypes.empty() or std::find(particleTypes.begin(),
                                               particleTypes.end(), part.type) != particleTypes.end()) {
            m = 1.0 * minima;
            k = 1.0 * kconstants;
            grad[0] = - 2 * k[0] * (x[0] - m[0]);
            grad[1] = - 2 * k[1] * (x[1] - m[1]);
            grad[2] = - 2 * k[2] * (x[2] - m[2]);
            force += scalefactor * grad;
        }
        return {force, torque};
    }
}
