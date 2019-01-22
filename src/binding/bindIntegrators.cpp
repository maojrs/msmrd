#include "binding.hpp"
#include "integrators/overdampedLangevin.hpp"
#include "integrators/overdampedLangevinMarkovSwitch.hpp"

namespace msmrd {
    /*
     * pyBinders for the c++ integrators classes
     */
    void bindIntegrators(py::module &m) {
        py::class_<overdampedLangevin, integrator>(m, "overdampedLangevin", "overdamped Langevin integrator (timestep, seed, "
                                                                "particlesbodytype (point, rod, rigidbody, "
                                                                "pointmix, rodmix or rigidbodymix) )")
                .def(py::init<double &, long &, std::string &>())
                .def_property_readonly("clock", &overdampedLangevin::getClock)
                .def("setKbT", &overdampedLangevin::setKbT)
                .def("setBoundary", &overdampedLangevin::setBoundary)
                .def("setExternalPotential", &overdampedLangevin::setExternalPotential)
                .def("setPairPotential", &overdampedLangevin::setPairPotential)
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
        py::bind_vector<std::vector<particle>>(m, "particleList", py::module_local());
        py::bind_vector<std::vector<particleMS>>(m, "particleMSList", py::module_local());
    }

}
