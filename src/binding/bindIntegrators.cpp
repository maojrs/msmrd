#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/pybind11.h>
#include "integrators/odLangevin.hpp"
#include "integrators/odLangevinMarkovSwitch.hpp"


// Needed to connect lists/arrays of particles in python with cpp integrator methods
PYBIND11_MAKE_OPAQUE(std::vector<particle>);
PYBIND11_MAKE_OPAQUE(std::vector<particleMS>);

namespace py = pybind11;

/*
 * pyBinders for the c++ integrators
 */
void bindIntegrators(py::module& m) {
    py::class_<odLangevin>(m, "odLangevin")
            .def(py::init<double&, long&, bool&>())
            .def_property_readonly("clock", &odLangevin::getClock)
            .def("setExternalPotential", &odLangevin::setExternalPotential)
            .def("setExternalRodPotential", &odLangevin::setExternalRodPotential)
                    //.def("evalExternalPotential", &odLangevin::evalExternalPotential)
            //.def("evalExternalForce", &odLangevin::evalExternalForce)
            .def("integrate", &odLangevin::integrate);

    py::class_<odLangevinMarkovSwitch<ctmsm>>(m, "odLangevinMarkovSwitch")
            .def(py::init<ctmsm&, double&, long&, bool&>())
            .def_property_readonly("clock", &odLangevinMarkovSwitch<ctmsm>::getClock)
            .def("setExternalPotential", &odLangevinMarkovSwitch<ctmsm>::setExternalPotential)
            .def("setExternalRodPotential", &odLangevinMarkovSwitch<ctmsm>::setExternalRodPotential)
                    //.def("evaluateExternalPotential", &odLangevin::evalExternalPotential)
            //.def("evalExternalForce", &odLangevin::evalExternalForce)
            .def("integrate", &odLangevinMarkovSwitch<ctmsm>::integrate);

    // Created c++ compatible particle list/vector/array of particles in python
    py::bind_vector<std::vector<particle>>(m, "particleList", py::module_local(false));
    py::bind_vector<std::vector<particleMS>>(m, "particleMSList", py::module_local(false));
}
