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
     */
    class patchyProteinMarkovSwitch : public patchyProtein {
    protected:
        double angularStrength = 2.0;
    public:
        /**
         * @param angularStrength give the angular strength of angular dependence of torque.
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
        double evaluate(const particle &part1, const particle &part2) override;

        std::array<vec3<double>, 4> forceTorque(const particle &part1, const particle &part2) override;


        // Additional auxiliary functions

        void setPotentialParameters() override;

        std::vector<std::vector<double>> forceTorquePyBind(particle &part1, particle &part2);

        std::tuple<vec3<double>, vec3<double>> calculatePlanes(const particle &part1, const particle &part2,
                                                               const std::vector<vec3<double>> patches1,
                                                               const std::vector<vec3<double>> patches2);
    };



}
