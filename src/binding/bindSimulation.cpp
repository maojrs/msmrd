#include "binding.hpp"
#include "particle.hpp"
#include "simulation.hpp"

namespace msmrd {
    /*
     * pyBinders for the c++ integrators classes
     */

void bindSimulation(py::module &m) {
        py::class_<simulation>(m, "simulation")
                .def(py::init<integrator &>())
                .def("run", &simulation::run);
        }
}