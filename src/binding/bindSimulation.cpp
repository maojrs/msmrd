#include "binding.hpp"
#include "simulation.hpp"

// Needed to connect lists/arrays of particles in python with cpp integrator methods
using particle = msmrd::particle;
using particleMS = msmrd::particleMS;
PYBIND11_MAKE_OPAQUE(std::vector<particle>);
PYBIND11_MAKE_OPAQUE(std::vector<particleMS>);

namespace msmrd {
    /*
     * pyBinders for the c++ integrators classes
     */
    void bindSimulation(py::module &m) {
        py::class_<simulation>(m, "simulation")
                .def(py::init<integrator &, std::vector<particle> &>())
                .def("run", &simulation::run);
        // Created c++ compatible particle list/vector/array of particles in python
        py::bind_vector<std::vector<particle>>(m, "particleList", py::module_local(false));
        py::bind_vector<std::vector<particleMS>>(m, "particleMSList", py::module_local(false));
    }
}