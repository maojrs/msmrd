//
// Created by maojrs on 9/7/18.
//
#include "binding.hpp"
#include "boundaries/boundary.hpp"
#include "boundaries/box.hpp"
#include "boundaries/sphere.hpp"



namespace msmrd {
    /*
     * pyBinders for the c++ boundaries classes
     */
    void bindBoundaries(py::module &m) {
        pybind11::class_<boundary>(m, "boundary");

        py::class_<box, boundary>(m, "box")
                .def(py::init<double &, double &, double &, std::string>())
                .def_property_readonly("boxsize", [](const box &currentbox) {
                    return vec2numpy(3, currentbox.boxsize);
                })
                .def_property_readonly("boundaryType", &box::getBoundaryType);

        py::class_<sphere, boundary>(m, "sphere")
                .def(py::init<double &, std::string>())
                .def_property_readonly("maxRadius", &sphere::getRadius)
                .def_property_readonly("boundaryType", &sphere::getBoundaryType);
    }

}
