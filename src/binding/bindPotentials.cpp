#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/cast.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>
#include <iostream>
#include "potentials/gaussians3D.hpp"
#include "potentials/harmonicRepulsion.hpp"

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"

namespace py = pybind11;

PYBIND11_MODULE(potentials, m) {

py::class_<gaussians3D>(m, "gaussians3D")
    .def(py::init<int&, double&, double&, long&>())
    .def("evaluate", &gaussians3D::evaluatePyBind)
    .def("force", &gaussians3D::forcePyBind);


#ifdef VERSION_INFO
m.attr("__version__") = VERSION_INFO;
#else
m.attr("__version__") = "dev";
#endif
}
#pragma clang diagnostic pop