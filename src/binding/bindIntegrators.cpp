#include "binding.hpp"
#include "integrators/overdampedLangevin.hpp"
#include "integrators/overdampedLangevinMarkovSwitch.hpp"
#include "integrators/msmrdIntegrator.hpp"

namespace msmrd {
    // Aliases for classes with long names.
    using msm = msmrd::discreteTimeMarkovStateModel;
    using ctmsm = msmrd::continuousTimeMarkovStateModel;
    using msmrdMSM = msmrd::msmrdMarkovStateModel;
    /*
     * pyBinders for the c++ integrators classes
     */
    void bindIntegrators(py::module &m) {
        py::class_<overdampedLangevin, integrator>(m, "overdampedLangevin", "overdamped Langevin integrator (timestep, seed, "
                                                                "particlesbodytype (point, rod, rigidbody, "
                                                                "pointmix, rodmix or rigidbodymix) )")
                .def(py::init<double &, long &, std::string &>())
                .def_property_readonly("clock", &overdampedLangevin::getClock)
                .def("setClock", &overdampedLangevin::setClock)
                .def("resetClock", &overdampedLangevin::resetClock)
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

        py::class_<msmrdIntegrator<ctmsm>, overdampedLangevinMarkovSwitch<ctmsm>>(m, "msmrdIntegrator",
                                                                              "Integrator for msmrd algorithm "
                                                                              "(timestep, seed, particlesbodytype, "
                                                                              "numParticleTypes, relativeDistanceCutOff,"
                                                                              "MSMlist, mainMarkovModel, "
                                                                              "positionOrientationPartition)")
                .def(py::init<double &, long &, std::string &, int &, double &, ctmsm &, msmrdMSM &>())
                .def("getRateFromKey", &msmrdIntegrator<ctmsm>::getRateFromKey)
                .def("setDiscretization", &msmrdIntegrator<ctmsm>::setDiscretization)
                .def("integrate", &msmrdIntegrator<ctmsm>::integrate);

        // Created c++ compatible particle list/vector/array of particles in python
        py::bind_vector<std::vector<particle>>(m, "particleList", py::module_local());
        py::bind_vector<std::vector<particleMS>>(m, "particleMSList", py::module_local());
    }

}
