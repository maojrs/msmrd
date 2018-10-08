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
        pybind11::class_<externalPotential<> >(m, "externalPotential");
        pybind11::class_<externalPotential<int> >(m, "externalMixPotential");
        pybind11::class_<externalPotential<vec3<double>>>(m, "externalRodPotential");
        pybind11::class_<externalPotential<vec3<double>, int>>(m, "externalRodMixPotential");
        pybind11::class_<externalPotential<quaternion<double>>>(m, "externalRigidBodyPotential");
        pybind11::class_<externalPotential<quaternion<double>, int>>(m, "externalRigidBodyMixPotential");

        /* Bind pair potential parent classes as instantiated templates
         * so inheritance works correctly. */
        pybind11::class_<pairPotential<> >(m, "pairPotential");
        pybind11::class_<pairPotential<int, int> >(m, "pairMixPotential");
        pybind11::class_<pairPotential<vec3<double>, vec3<double>>>(m, "pairRodPotential");
        pybind11::class_<pairPotential<vec3<double>, vec3<double>, int, int>>(m, "pairRodMixPotential");
        pybind11::class_<pairPotential<quaternion<double>, quaternion<double>>>(m, "pairRigidBodyPotential");
        pybind11::class_<pairPotential<quaternion<double>, quaternion<double>, int, int>>
                (m, "pairRigidBodyMixPotential");
    }

}
