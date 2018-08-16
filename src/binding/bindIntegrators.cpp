#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/cast.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>
#include <iostream>
#include "integrators/odLangevin.hpp"
#include "integrators/odLangevinMarkovSwitch.hpp"

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"

// Needed to connect lists/arrays of particles in python with cpp integrator methods
PYBIND11_MAKE_OPAQUE(std::vector<particle>);
PYBIND11_MAKE_OPAQUE(std::vector<particleMS>);

namespace py = pybind11;

PYBIND11_MODULE(integrators, m) {


/*
 * pyBinders for all the relevant c++ classes
 */
py::class_<odLangevin>(m, "odLangevin")
    .def(py::init<double&, long&, bool&>())
    .def_property("clock", &odLangevin::getClock, nullptr)
    .def("setPotential", &odLangevin::setPotential)
    .def("integrate", &odLangevin::integrate)
    .def("integrateList", &odLangevin::integrateList);

py::class_<odLangevinMarkovSwitch<ctmsm>>(m, "odLangevinMarkovSwitch")
    .def(py::init<ctmsm&, double&, long&, bool&>())
    .def_property("clock", &odLangevinMarkovSwitch<ctmsm>::getClock, nullptr)
    .def("integrate", &odLangevinMarkovSwitch<ctmsm>::integrate)
    .def("integrateList", &odLangevinMarkovSwitch<ctmsm>::integrateList);

// Created c++ compatible particle list/vector/array of particles in python
py::bind_vector<std::vector<particle>>(m, "particleList");
py::bind_vector<std::vector<particleMS>>(m, "particleMSList");

#ifdef VERSION_INFO
m.attr("__version__") = VERSION_INFO;
#else
m.attr("__version__") = "dev";
#endif
}
#pragma clang diagnostic pop