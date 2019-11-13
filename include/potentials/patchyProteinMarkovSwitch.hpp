//
// Created by maojrs on 10/23/19.
//

#pragma once
#include "potentials/patchyProtein.hpp"

namespace msmrd {
    /*
     * Declaration of potential function for patchy protein model with angular dependence. This class is a
     * child from patchyProtein class; it accepts quaternion-based orientation (rigidBody particle type). The main
     * idea is to have a potential with several different patches (patchy protein) that have a fixed angular
     * binding. It further allows that one of the patches can be turned on and off depending on the MSM state.
     *
     * Warning: since these classes are highly dependent on a given implementation, they have many hard coded
     * values. If a new implemantation is desired, we recommend to create a new child class of patchyProtein
     * in this file using the classes here as reference (see setPotentialParameters() and calculatePlanes()).
     */
    class patchyProteinMarkovSwitch : public patchyProtein {
    protected:
        double angularStrength = 2.0;
        double minimumR = 1.25;
    public:
        /**
         * @param angularStrength give the angular strength of angular dependence of torque.
         * @param minimumR minimum distance that the particles need to be away, so normal MSM behavior is active.
         */

        // Inherit parent class contructor
        using patchyProtein::patchyProtein;

        // Additional constructors (in case angular strength is determined)
        patchyProteinMarkovSwitch(double sigma, double strength, double angularStrength,
                      std::vector<vec3<double>> patchesCoordinatesA,
                      std::vector<vec3<double>> patchesCoordinatesB);

        patchyProteinMarkovSwitch(double sigma, double strength, double angularStrength,
                      std::vector<std::vector<double>> patchesCoordinatesA,
                      std::vector<std::vector<double>> patchesCoordinatesB);

        /* Note evaluate and forceTorque functions do not override the ones of patchyProtein since these
         * ones take particle as arguments instead of particle. Therefore one must be careful the integrator
         * uses particle if we want these functions to be used. */
        double evaluate(particle &part1, particle &part2) override;

        std::array<vec3<double>, 4> forceTorque(particle &part1, particle &part2) override;


        // Additional auxiliary functions

        void setPotentialParameters() override;

        void enableDisableMSM(vec3<double>relPosition, particle &part1, particle &part2);

        std::tuple<vec3<double>, vec3<double>> calculatePlanes(particle &part1, particle &part2,
                                                               const std::vector<vec3<double>> patches1,
                                                               const std::vector<vec3<double>> patches2);
    };



}
