//
// Created by maojrs on 3/18/21.
//
#pragma once
#include "potentials/patchyProtein.hpp"

namespace msmrd{

    /*
     * Declaration of potential function for patchy protein model of the MAPK pathway. This class is a child from the
     * general pair potentials class that accepts quaternion-based orientation (rigidBody particle type).
     * Some of its functionality is based on the patchy particles class. The main difference is that this class
     * can take different type of patches. Note the potential will depend on the position of both particles,
     * their orientation (quaternion<double>, orientvector<double>), and their types (int, int). Note the orientation
     * of the ligand (kinase or phosphotase) only depends on the orientation vector (has one free orientational degree
     * of freedom).
     *
     * WORKING TO MODIFY THIS POTENTIAL
     */
    class patchyProteinMAPK : public patchyProtein {
    private:
        // Need to call this function from constructor to define parameters
        void setPotentialParameters();

    public:

        using patchyProtein::patchyProtein;

    };

}