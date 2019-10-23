//
// Created by maojrs on 10/23/19.
//

#pragma once
#include "potentials/patchyProtein.hpp"

namespace msmrd{
    /*
     * Declaration of potential function for patchy protein model with angular dependence. This class is a
     * child from patchyProtein class; it accepts quaternion-based orientation (rigidBody particle type). The main
     * idea is to have a potential with several different patches (patchy protein) that have a fixed angular
     * binding. It further allows that one of the patches can be turned on and off depending on the MSM state.
     */
    class patchyProteinAngular : public patchyProtein {

        // Inherit parent class contructor
        using patchyProtein::patchyProtein;

        double evaluate(const particle &part1, const particle &part2) override;

        std::array<vec3<double>, 4>
        forceTorque(const particle &part1, const particle &part2) override;

    };
}
