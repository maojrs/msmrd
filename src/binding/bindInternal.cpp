//
// Created by maojrs on 10/8/18.
//
#include "binding.hpp"
#include "potentials/potentials.hpp"
#include "boundaries/boundary.hpp"



namespace msmrd {
    /*
     * pyBinders for the c++ parent classes that should remain hidden under _internal in python
     */
    void bindInternal(py::module &m) {

        /* Bind boundaries parent class*/
        pybind11::class_<boundary>(m, "boundary");

        /* Bind external potential parent classes as instantiated templates
         * so inheritance works correctly. */
        pybind11::class_<externalPotential >(m, "externalPotential");

        /* Bind pair potential parent classes as instantiated templates
         * so inheritance works correctly. */
        pybind11::class_<pairPotential >(m, "pairPotential");
                (m, "pairRigidBodyMixPotential");
    }

}
