#include "binding.hpp"
#include "integrators/integratorMAPK.hpp"
#include "integrators/integratorMoriZwanzig.hpp"
#include "integrators/langevin.hpp"
#include "integrators/overdampedLangevin.hpp"
#include "integrators/overdampedLangevinMarkovSwitch.hpp"
#include "integrators/overdampedLangevinSelective.hpp"
#include "integrators/msmrdIntegrator.hpp"
#include "integrators/msmrdPatchyProtein.hpp"
#include "integrators/msmrdMultiParticleIntegrator.hpp"
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

        py::class_<langevin, integrator>(m, "langevin", "Langevin integrator (timestep, seed, "
                                                        "particlesbodytype (point, rod, rigidbody, "
                                                        "pointmix, rodmix or rigidbodymix), frictionCoefficient,"
                                                        "integratorScheme (BAOAB default) )")
                .def(py::init<double &, long &, std::string &, double &>())
                .def(py::init<double &, long &, std::string &, double &, std::string &>())
                .def("integrate", &langevin::integrate);

        py::class_<overdampedLangevinSelective, overdampedLangevin>(m, "overdampedLangevinSelective", "overdamped "
                                                                    "Langevin integrator with selective active patches."
                                                                    " Special for multiparticle simulations w/patchy "
                                                                    "particles that avoid triple bindings. (timestep, "
                                                                    "seed, particlesbodytype (point, rod, rigidbody, "
                                                                            "pointmix, rodmix or rigidbodymix) )")
                .def(py::init<double &, long &, std::string &>())
                .def("integrate", &overdampedLangevinSelective::integrate)
                .def("updateParticleCompounds", &overdampedLangevinSelective::updateParticleCompounds)
                .def("findClosedBindingLoops", &overdampedLangevinSelective::findClosedBindingLoops)
                .def("getCompoundSize", &overdampedLangevinSelective::getCompoundSize);

        py::class_<integratorMoriZwanzig, langevin>(m, "integratorMoriZwanzig", "Specialized Langevin integrator "
                                                        "for MoriZwanzig application (timestep, seed, "
                                                        "particlesbodytype (point, rod, rigidbody, "
                                                        "pointmix, rodmix or rigidbodymix) )")
                .def(py::init<double &, long &, std::string &, double &>())
                .def("integrate", &integratorMoriZwanzig::integrate)
                .def("setDistinguishedTypes", &integratorMoriZwanzig::setDistinguishedTypes);

        py::class_<integratorMAPK, overdampedLangevin>(m, "integratorMAPK", "integrator for MAPK (timestep, seed, "
                                                                            "particlesbodytype (point, rod, rigidbody, "
                                                                            "pointmix, rodmix or rigidbodymix), "
                                                                            "anglePatches, mapkIndex, kinaseIndex, "
                                                                            "phosIndex )")
                .def(py::init<double &, long &, std::string &, double &, double &, double &, double &,
                        std::vector<int> &, std::vector<int> &, std::vector<int> &>())
                .def("integrate", &integratorMAPK::integrate)
                .def("disableMAPK", &integratorMAPK::disableMAPK);

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
                .def("setRecordEventLog", &msmrdIntegrator<ctmsm>::setRecordEventLog)
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


        py::class_<msmrdPatchyProtein2, msmrdPatchyProtein>(m, "msmrdPatchyProtein2",
                                                               "Integrator for msmrd algorithm, specialized for"
                                                               "patchyProtein2 implementation (timestep, seed, "
                                                               "particlesbodytype, numParticleTypes, "
                                                               "relativeDistanceCutOff,"
                                                               "MSMlist, mainMarkovModel, "
                                                               "positionOrientationPartition)")
                .def(py::init<double &, long &, std::string &, int &, std::array<double,2> &, ctmsm &, msmrdMSM &>())
                .def(py::init<double &, long &, std::string &, int &, std::array<double,2> &, std::vector<ctmsm> &,
                        msmrdMSM &>());

        py::class_<msmrdMultiParticleIntegrator<ctmsm>, msmrdIntegrator<ctmsm>>(m, "msmrdMultiParticleIntegrator",
                                                               "Integrator for msmrd algorithm, specialized for"
                                                               "multiparticle intergation of dimer molecules"
                                                               "(timestep, seed, "
                                                               "particlesbodytype, numParticleTypes, "
                                                               "radialBounds, MSMlist, mainMarkovModel),"
                                                               "DlistCompound, DrotlistCompound ")
                .def(py::init<double &, long &, std::string &, int &, std::array<double,2> &, ctmsm &, msmrdMSM &,
                        std::vector<double> &, std::vector<double> &>())
                .def(py::init<double &, long &, std::string &, int &, std::array<double,2> &, std::vector<ctmsm> &,
                        msmrdMSM, std::vector<double> &, std::vector<double> &>())
                .def("findClosedBindingLoops", &msmrdMultiParticleIntegrator<ctmsm>::findClosedBindingLoops, "finds "
                                      "closed binding loops in all compounds")
                .def("getNumberOfBindingsInCompound", &msmrdMultiParticleIntegrator<ctmsm>::getNumberOfBindingsInCompound,
                     "gets number of bindings in a give compound")
                .def("getBindingsInCompound", &msmrdMultiParticleIntegrator<ctmsm>::getBindingsInCompound,
                     "gets bindings in a give compound")
                .def("getCompoundSize", &msmrdMultiParticleIntegrator<ctmsm>::getCompoundSize,
                     "gets compound size");


        // Created c++ compatible particle list/vector/array of particles in python
        py::bind_vector<std::vector<particle>>(m, "particleList", py::module_local());

    }

}
