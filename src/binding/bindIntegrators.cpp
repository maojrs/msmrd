#include "binding.hpp"
#include "integrators/odLangevin.hpp"
#include "integrators/odLangevinMarkovSwitch.hpp"


// Needed to connect lists/arrays of particles in python with cpp integrator methods
using particle = msmrd::particle;
using particleMS = msmrd::particleMS;
PYBIND11_MAKE_OPAQUE(std::vector<particle>);
PYBIND11_MAKE_OPAQUE(std::vector<particleMS>);

namespace msmrd {
    /*
     * pyBinders for the c++ integrators classes
     */
    void bindIntegrators(py::module &m) {
        py::class_<odLangevin>(m, "odLangevin")
                .def(py::init<double &, long &, bool &>())
                .def_property_readonly("clock", &odLangevin::getClock)
                .def("setKbT", &odLangevin::setKbT)
                .def("setBoundary", &odLangevin::setBoundary)
                .def("setExternalPotential", &odLangevin::setExternalPotential)
                .def("setExternalRodPotential", &odLangevin::setExternalRodPotential)
                .def("setPairPotential", &odLangevin::setPairPotential)
                .def("setPairRodPotential", &odLangevin::setPairRodPotential)
                        //.def("evalExternalPotential", &odLangevin::evalExternalPotential)
                        //.def("evalExternalForce", &odLangevin::evalExternalForce)
                .def("integrate", &odLangevin::integrate);

        py::class_<odLangevinMarkovSwitch<ctmsm>, odLangevin>(m, "odLangevinMarkovSwitch")
                .def(py::init<ctmsm &, double &, long &, bool &>())
                .def("integrate", &odLangevinMarkovSwitch<ctmsm>::integrate);

        // Created c++ compatible particle list/vector/array of particles in python
        py::bind_vector<std::vector<particle>>(m, "particleList", py::module_local(false));
        py::bind_vector<std::vector<particleMS>>(m, "particleMSList", py::module_local(false));
    }

}
