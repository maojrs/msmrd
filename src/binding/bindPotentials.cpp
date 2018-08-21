#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/cast.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>
#include "potentials/potentials.hpp"
#include "potentials/gaussians3D.hpp"
#include "potentials/harmonicRepulsion.hpp"


namespace py = pybind11;

/*
 * pyBinders for the c++ potentials
 */
void bindPotentials(py::module& m) {
    pybind11::class_<externalPotential>(m, "externalPotential");

    py::class_<gaussians3D, externalPotential>(m, "gaussians3D")
        .def(py::init<int&, double&, double&, long&>())
        .def("evaluate", &gaussians3D::evaluatePyBind)
        .def("force", &gaussians3D::forcePyBind);
}
