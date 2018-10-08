#include "binding.hpp"
#include "integrators/overdampedLangevin.hpp"
#include "integrators/overdampedLangevinMarkovSwitch.hpp"


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
        py::class_<overdampedLangevin>(m, "overdampedLangevin", "overdamped Langevin integrator (timestep, seed, "
                                                                "particlesbodytype (point, rod, rigidbody, "
                                                                "pointmix, rodmix or rigidbodymix) )")
                .def(py::init<double &, long &, std::string &>())
                .def_property_readonly("clock", &overdampedLangevin::getClock)
                .def("setKbT", &overdampedLangevin::setKbT)
                .def("setBoundary", &overdampedLangevin::setBoundary)
                .def("setExternalPotential", &overdampedLangevin::setExternalPotential)
                .def("setExternalMixPotential", &overdampedLangevin::setExternalMixPotential)
                .def("setExternalRodPotential", &overdampedLangevin::setExternalRodPotential)
                .def("setExternalRodMixPotential", &overdampedLangevin::setExternalRodMixPotential)
                .def("setExternalRigidBodyPotential", &overdampedLangevin::setExternalRigidBodyPotential)
                .def("setExternalRigidBodyMixPotential", &overdampedLangevin::setExternalRigidBodyMixPotential)
                .def("setPairPotential", &overdampedLangevin::setPairPotential)
                .def("setPairMixPotential", &overdampedLangevin::setPairMixPotential)
                .def("setPairRodPotential", &overdampedLangevin::setPairRodPotential)
                .def("setPairRodMixPotential", &overdampedLangevin::setPairRodMixPotential)
                .def("setPairRigidBodyPotential", &overdampedLangevin::setPairRigidBodyPotential)
                .def("setPairRigidBodyMixPotential", &overdampedLangevin::setPairRigidBodyMixPotential)
                .def("integrate", &overdampedLangevin::integrate);

        py::class_<overdampedLangevinMarkovSwitch<ctmsm>, overdampedLangevin>(m, "overdampedLangevinMarkovSwitch",
                                                                              "overdamped Langevin integrator with "
                                                                              "Markov switch (Markov model, timestep,"
                                                                              " seed, particlesbodytype (point, "
                                                                              "rod, rigidbody, pointmix, rodmix or "
                                                                              "rigidbodymix) )")
                .def(py::init<ctmsm &, double &, long &, std::string &>())
                .def("integrate", &overdampedLangevinMarkovSwitch<ctmsm>::integrate);

        // Created c++ compatible particle list/vector/array of particles in python
        py::bind_vector<std::vector<particle>>(m, "particleList", py::module_local(false));
        py::bind_vector<std::vector<particleMS>>(m, "particleMSList", py::module_local(false));
    }

}
