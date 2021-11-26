#include "binding.hpp"
#include "simulation.hpp"

namespace msmrd {
    /*
     * pyBinders for the c++ integrators classes
     */

void bindSimulation(py::module &m) {
        py::class_<simulation>(m, "simulation")
                .def(py::init<integrator &>())
                .def("setEquilibrationSteps", &simulation::setEquilibrationSteps)
                .def("setOutputEnergyTemperature", &simulation::setOutputEnergyTemperature)
                .def("setDistinguishedTypes", &simulation::setDistinguishedTypes)
                .def("run", &simulation::run);
        }
}
