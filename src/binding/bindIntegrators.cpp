#include "binding.hpp"
#include "integrators/overdampedLangevin.hpp"
#include "integrators/overdampedLangevinMarkovSwitch.hpp"
#include "integrators/msmrdIntegrator.hpp"
#include "integrators/msmrdPatchyProtein.hpp"
#include "discretizations/positionOrientationPartition.hpp"


namespace msmrd {
    // Aliases for classes with long names.
    using msm = msmrd::discreteTimeMarkovStateModel;
    using ctmsm = msmrd::continuousTimeMarkovStateModel;
    using msmrdMSM= msmrd::msmrdMarkovModel;


    /*
     * pyBinders for the c++ integrators classes
     */
    void bindIntegrators(py::module &m) {
        py::class_<overdampedLangevin, integrator>(m, "overdampedLangevin", "overdamped Langevin integrator (timestep, seed, "
                                                                "particlesbodytype (point, rod, rigidbody, "
                                                                "pointmix, rodmix or rigidbodymix) )")
                .def(py::init<double &, long &, std::string &>())
                .def("integrate", &overdampedLangevin::integrate);

        py::class_<overdampedLangevinMarkovSwitch<ctmsm>, overdampedLangevin>(m, "overdampedLangevinMarkovSwitch",
                                                                              "overdamped Langevin integrator with "
                                                                              "Markov switch (Markov model, timestep,"
                                                                              " seed, particlesbodytype (point, "
                                                                              "rod, rigidbody, pointmix, rodmix or "
                                                                              "rigidbodymix) )")
                .def(py::init<ctmsm &, double &, long &, std::string &>())
                .def(py::init<std::vector<ctmsm> &, double &, long &, std::string &>())
                .def("integrate", &overdampedLangevinMarkovSwitch<ctmsm>::integrate);

        /* Note this binding uses method overloading, and it therefore requires to explicitly state
         * the input arguments when functions are overloaded. */
        py::class_<msmrdIntegrator<ctmsm>, overdampedLangevinMarkovSwitch<ctmsm>>(m, "msmrdIntegrator",
                                                                                  "Integrator for msmrd algorithm "
                                                                                  "(timestep, seed, particlesbodytype, "
                                                                                  "numParticleTypes, "
                                                                                  "relativeDistanceCutOff,"
                                                                                  "MSMlist, mainMarkovModel, "
                                                                                  "positionOrientationPartition)")
                .def(py::init<double &, long &, std::string &, int &, std::array<double,2> &, ctmsm &, msmrdMSM &>())
                .def(py::init<double &, long &, std::string &, int &, std::array<double,2> &, std::vector<ctmsm> &,
                        msmrdMSM &>())
                .def("setDiscretization", (void (msmrdIntegrator<ctmsm>::*)
                        (std::shared_ptr<spherePartition>&)) &msmrdIntegrator<ctmsm>::setDiscretization,
                        "Sets position only discretization")
                .def("setDiscretization", (void (msmrdIntegrator<ctmsm>::*)
                        (std::shared_ptr<positionOrientationPartition>&)) &msmrdIntegrator<ctmsm>::setDiscretization,
                        "Sets full position orientation discretization")
                .def("printEventLog", &msmrdIntegrator<ctmsm>::printEventLog)
                .def("integrate", &msmrdIntegrator<ctmsm>::integrate);


        py::class_<msmrdPatchyProtein, msmrdIntegrator<ctmsm>>(m, "msmrdPatchyProtein",
                                                                  "Integrator for msmrd algorithm, specialized for"
                                                                  "patchyProtein implementation (timestep, seed, "
                                                                  "particlesbodytype, numParticleTypes, "
                                                                  "relativeDistanceCutOff,"
                                                                  "MSMlist, mainMarkovModel, "
                                                                  "positionOrientationPartition)")
                .def(py::init<double &, long &, std::string &, int &, std::array<double,2> &, ctmsm &, msmrdMSM &>())
                .def(py::init<double &, long &, std::string &, int &, std::array<double,2> &, std::vector<ctmsm> &,
                        msmrdMSM &>());


        // Created c++ compatible particle list/vector/array of particles in python
        py::bind_vector<std::vector<particle>>(m, "particleList", py::module_local());

    }

}
