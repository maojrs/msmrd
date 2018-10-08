//
// Created by maojrs on 9/7/18.
//

#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;

namespace msmrd {
    /*
     * Define functions to bind modules and submodules (implemented in src/binding/bind****.py)
     */
    void bindBoundaries(py::module&);
    void bindIntegrators(py::module&);
    void bindMarkovModels(py::module&);
    void bindParticles(py::module&);
    void bindPotentials(py::module&);
    void bindTrajectories(py::module&);

    // Function template to transform vectors to numpy for pyBindings
    template<typename ndvec>
    py::array_t<double> vec2numpy(int size, ndvec v) {
        auto result = py::array_t<double>(size);
        py::buffer_info buf = result.request();
        auto *ptr = (double *) buf.ptr;
        for (size_t idx = 0; idx < size; idx++)
            ptr[idx] = v[idx];
        return result;
    }
}
