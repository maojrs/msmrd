//
// Created by maojrs on 2/2/22.
//

#include "particle.hpp"
#include "vec3.hpp"
#include "potentials/bistable.hpp"

namespace msmrd {
    /**
     * Class implementation for harmonic external potential
     */

    /* Alternative constructor to specify location of minimas and standard deviation explicitly. */
    bistable::bistable(double minimaDist, std::vector<double> kconstants_in, double scalefactor)
            : minimaDist(minimaDist), scalefactor(scalefactor) {
        if (kconstants_in.size() !=3 ){
            throw std::invalid_argument("Kconstants should be three-dimensional arrays");
        }
        kconstants = vec3<double>{kconstants_in[0], kconstants_in[1], kconstants_in[2]};
    }

    /* Alternative constructor to specify location of minimas and standard deviation explicitly, as well as
     * the specific particle types on which the potential will act on. */
    bistable::bistable(double minimaDist, std::vector<double> kconstants_in,
                       std::vector<int> partTypes, double scalefactor) :
            bistable(minimaDist, kconstants_in, scalefactor)  {
        particleTypes = partTypes;
    }


    // Returns value of potential at position x
    double bistable::evaluate(particle &part) {
        vec3<double> x = part.position;
        vec3<double> k;
        double bistablePot;
        if (particleTypes.empty() or std::find(particleTypes.begin(),
                                               particleTypes.end(), part.type) != particleTypes.end()) {
            k = 1.0 * kconstants;
            bistablePot = k[0] * std::pow(1 - std::pow(x[0]/minimaDist,2),2) +
                          k[1] * std::pow(x[1], 2) +
                          k[2] * std::pow(x[2], 2) ;
        }
        return scalefactor * bistablePot;
    };


    // Returns minus gradient of potential (force) at position x and zero torque
    std::array<vec3<double>, 2> bistable::forceTorque(particle &part) {
        vec3<double> x = part.position;
        vec3<double> force = vec3<double>(0, 0, 0);
        vec3<double> torque = vec3<double>(0, 0, 0);
        vec3<double> k, grad;
        if (particleTypes.empty() or std::find(particleTypes.begin(),
                                               particleTypes.end(), part.type) != particleTypes.end()) {
            k = 1.0 * kconstants;
            grad[0] = - k[0] * 4 * x[0] * (std::pow(x[0],2) - std::pow(minimaDist, 2))/std::pow(minimaDist,4);
            grad[1] = - 2 * k[1] * (x[1]);
            grad[2] = - 2 * k[2] * (x[2]);
            force += scalefactor * grad;
        }
        return {force, torque};
    }


    /*
     * Implementations of bistable2, an alternative version
     */

    /* Alternative constructor to specify location of minimas and standard deviation explicitly. */
    bistable2::bistable2(std::vector<double> parameters_in, double scalefactor)
            : scalefactor(scalefactor) {
        if (parameters_in.size() !=4 ){
            throw std::invalid_argument("Parameters should be a four-dimensional arrays");
        }
        parameters = parameters_in;
    }

    /* Alternative constructor to specify location of minimas and standard deviation explicitly, as well as
     * the specific particle types on which the potential will act on. */
    bistable2::bistable2(std::vector<double> parameters_in,
                       std::vector<int> partTypes, double scalefactor) :
            bistable2(parameters_in, scalefactor)  {
        particleTypes = partTypes;
    }


    // Returns value of potential at position x
    double bistable2::evaluate(particle &part) {
        vec3<double> x = part.position;
        double bistablePot;
        std::vector<double>& p = parameters;
        if (particleTypes.empty() or std::find(particleTypes.begin(),
                                               particleTypes.end(), part.type) != particleTypes.end()) {
            bistablePot = p[0] * std::pow(x[0],4) - p[1] * std::pow(x[0],2) +
                          p[2] * std::pow(x[1], 2) + p[3] * std::pow(x[2], 2);
        }
        return scalefactor * bistablePot;
    };


    // Returns minus gradient of potential (force) at position x and zero torque
    std::array<vec3<double>, 2> bistable2::forceTorque(particle &part) {
        vec3<double> x = part.position;
        vec3<double> force = vec3<double>(0, 0, 0);
        vec3<double> torque = vec3<double>(0, 0, 0);
        vec3<double> grad;
        std::vector<double>& p = parameters;
        if (particleTypes.empty() or std::find(particleTypes.begin(),
                                               particleTypes.end(), part.type) != particleTypes.end()) {
            grad[0] = - 4 * p[0] * std::pow(x[0],3) + 2 * p[1] * x[0];
            grad[1] = - 2 * p[2] * (x[1]);
            grad[2] = - 2 * p[3] * (x[2]);
            force += scalefactor * grad;
        }
        return {force, torque};
    }

}
