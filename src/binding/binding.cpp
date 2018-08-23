#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/cast.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>
#include <numeric>
#include <iostream>
#include <random>
#include "msm.hpp"
#include "vec3.hpp"
#include "quaternion.hpp"
#include "particle.hpp"
#include "simulation.hpp"


#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"


namespace py = pybind11;

// Define functions to bind submodules (implemented in src/binding/bind****.py)
void bindIntegrators(py::module&);
void bindPotentials(py::module&);


/* Convert c++ vectors/arrays to numpy arrays for pybind */
template <typename ndvec>
py::array_t<double> vec2numpy(int size, ndvec v) {
    auto result = py::array_t<double>(size);
    py::buffer_info buf = result.request();
    auto *ptr = (double *) buf.ptr;
    for (size_t idx=0; idx<size; idx++)
        ptr[idx] = v[idx];
    return result;
}

namespace py = pybind11;

PYBIND11_MODULE(msmrd2binding, m) {

    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: msmrd2

        .. autosummary::
           :toctree: _generate

           particle
           msm
           ctmsm
    )pbdoc";

    /**
     * Call functions to bind msmrd2 submodules (functions defined above
     * and implemented in src/bind****.py)
     */
    {
        // Load integrators submodule
        auto integratorsSubmodule = m.def_submodule("integrators");
        bindIntegrators(integratorsSubmodule);
        // Load potentials submodule
        auto potentialsSubmodule = m.def_submodule("potentials");
        bindPotentials(potentialsSubmodule);
    }


    /*
     * pyBinders for the main c++ classes
     */
    py::class_<particle>(m, "particle")
            .def(py::init<double&, double&, std::string&, std::vector<double>&, std::vector<double>& >())
            .def_property_readonly("ID", &particle::getID)
            .def_property_readonly("D", &particle::getD)
            .def_property_readonly("Drot", &particle::getDrot)
            .def_property_readonly("bodytype", &particle::getBodyType)
            .def_property_readonly("position", [](const particle &part) {
                return vec2numpy(3,part.position);
            })
            .def_property_readonly("orientation", [](const particle &part) {
                return vec2numpy(4,part.orientation);
            })
            .def("setD", &particle::setD)
            .def("setDrot", &particle::setDrot)
            .def("setBodyType", &particle::setBodyType)
            .def("setPosition", &particle::setPositionPyBind)
            .def("setOrientationVec", &particle::setOrientationVecPyBind)
            .def("setOrientation", &particle::setOrientationPyBind);

    py::class_<particleMS, particle>(m, "particleMS")
            .def(py::init<int&, int&, double&, double&, std::string&, std::vector<double>&, std::vector<double>& >())
            .def_property_readonly("type", &particleMS::getType)
            .def_property_readonly("state", &particleMS::getState)
            .def_property_readonly("lagtime", &particleMS::getLagtime)
            .def("setState", &particleMS::setState)
            .def("setType", &particleMS::setType)
            .def("setLagtime", &particleMS::setLagtime);

    py::class_<msm>(m, "msm")
            .def(py::init<int&, std::vector<std::vector<double>>&, double&, long&>())
            .def_property_readonly("ID", &msm::getID)
            .def_property_readonly("nstates", &msm::getNstates)
            .def_property_readonly("lagtime", &msm::getLagtime)
            .def_property_readonly("tmatrix", &msm::getTmatrix)
            .def_property_readonly("D", [](const msm &currentmsm) {
                return vec2numpy(currentmsm.nstates,currentmsm.Dlist);
            })
            .def_property_readonly("Drot", [](const msm &currentmsm) {
                return vec2numpy(currentmsm.nstates,currentmsm.Drotlist);
            })
            .def("setD", &msm::setD)
            .def("setDrot", &msm::setDrot)
            .def("propagate", &msm::propagate);

    py::class_<ctmsm>(m, "ctmsm")
            .def(py::init<int&, std::vector<std::vector<double>>&, long&>())
            .def_property_readonly("ID", &ctmsm::getID)
            .def_property_readonly("nstates", &ctmsm::getNstates)
            .def_property_readonly("lagtime", &ctmsm::getLagtime)
            .def_property_readonly("tmatrix", &ctmsm::getTmatrix)
            .def_property_readonly("D", [](const ctmsm &currentmsm) {
                return vec2numpy(currentmsm.nstates,currentmsm.Dlist);
            })
            .def_property_readonly("Drot", [](const msm &currentmsm) {
                return vec2numpy(currentmsm.nstates,currentmsm.Drotlist);
            })
            .def("setD", &ctmsm::setD)
            .def("setDrot", &ctmsm::setDrot)
            .def("propagate", &ctmsm::propagate);


#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
#pragma clang diagnostic pop