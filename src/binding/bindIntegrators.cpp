#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/pybind11.h>
#include "integrators/odLangevin.hpp"
#include "integrators/odLangevinMarkovSwitch.hpp"

namespace py = pybind11;

/*
 * pyBinders for the c++ integrators
 */
void bindIntegrators(py::module& m) {
    py::class_<odLangevin>(m, "odLangevin")
            .def(py::init<double&, long&, bool&>())
            .def_property_readonly("clock", &odLangevin::getClock)
            .def("setExternalPotential", &odLangevin::setExternalPotential)
            //.def("evalExternalPotential", &odLangevin::evalExternalPotential)
            //.def("evalExternalForce", &odLangevin::evalExternalForce)
            .def("integrate", &odLangevin::integrate)
            .def("integrateList", &odLangevin::integrateList);

    py::class_<odLangevinMarkovSwitch<ctmsm>>(m, "odLangevinMarkovSwitch")
            .def(py::init<ctmsm&, double&, long&, bool&>())
            .def_property_readonly("clock", &odLangevinMarkovSwitch<ctmsm>::getClock)
            .def("setExternalPotential", &odLangevin::setExternalPotential)
            //.def("evaluateExternalPotential", &odLangevin::evalExternalPotential)
            //.def("evalExternalForce", &odLangevin::evalExternalForce)
            .def("integrate", &odLangevinMarkovSwitch<ctmsm>::integrate)
            .def("integrateList", &odLangevinMarkovSwitch<ctmsm>::integrateList);
}
