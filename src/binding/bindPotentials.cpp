#include <pybind11/pybind11.h>
#include "potentials/gaussians3D.hpp"
#include "potentials/harmonicRepulsion.hpp"


namespace py = pybind11;

/*
 * pyBinders for the c++ potentials
 */
void bindPotentials(py::module& m) {
    py::class_<gaussians3D>(m, "gaussians3D")
        .def(py::init<int&, double&, double&, long&>())
        .def("evaluate", &gaussians3D::evaluatePyBind)
        .def("force", &gaussians3D::forcePyBind);
}
