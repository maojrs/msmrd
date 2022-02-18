//
// Created by maojrs on 2/2/22.
//

#pragma once
#include "vec3.hpp"
#include "potentials/potentials.hpp"

namespace msmrd {
    /*
     * Bistable external potential in 3D corresponding to scaleFactor ( k1 (1-(x/minimaDist)^2)^2 + k2 y^2 + k3 z^2).
     */
    class bistable : public externalPotential {
    public:
        double scalefactor = 1;
        long seed = 0;
        std::vector<int> particleTypes;
        double minimaDist;
        vec3<double> kconstants;

        /**
         * @param scalefactor adjusts strength of the potential muliplying by a constant
         * @param seed variable for random number generation (Note seed = -1 corresponds to random device)
         * @param particleTypes if assigned values, the external potential only acts on the specified
         * particle types.
         * @param minimaDist distance of each of the two minimas to the origin
         * @param kconstants 3d vectors of k constant in each dimension
         */
        bistable(double minimaDist, std::vector<double> kconstants, double scalefactor);

        bistable(double minimaDist, std::vector<double> kconstants,
                 std::vector<int> particleTypes, double scalefactor);

        double evaluate(particle &part) override;

        std::array<vec3<double>, 2> forceTorque(particle &part) override;
    };


   /*
    * Alternative version of bistable external potential in 3D corresponding to:
    * scaleFactor scalefactor( ax^4 - bx^2 + cy^2 + dz^2).
    */
    class bistable2 : public externalPotential {
    public:
        double scalefactor = 1;
        long seed = 0;
        std::vector<int> particleTypes;
        std::vector<double> parameters;
        /**
         * @param scalefactor adjusts strength of the potential muliplying by a constant
         * @param seed variable for random number generation (Note seed = -1 corresponds to random device)
         * @param particleTypes if assigned values, the external potential only acts on the specified
         * particle types.
         * @param parameters , paramaters a,b,c,d to determine form of potential
         */
        bistable2(std::vector<double> parameters, double scalefactor);

        bistable2(std::vector<double> parameters, std::vector<int> particleTypes, double scalefactor);

        double evaluate(particle &part) override;

        std::array<vec3<double>, 2> forceTorque(particle &part) override;
    };

}

