#include "binding.hpp"
#include "potentials/potentials.hpp"
#include "potentials/bistable.hpp"
#include "potentials/combinedPotentials.hpp"
#include "potentials/dipole.hpp"
#include "potentials/gaussians3D.hpp"
#include "potentials/gayBerne.hpp"
#include "potentials/harmonic.hpp"
#include "potentials/harmonicRepulsion.hpp"
#include "potentials/lennardJones.hpp"
#include "potentials/pairBistable.hpp"
#include "potentials/patchyParticle.hpp"
#include "potentials/patchyParticleAngular.hpp"
#include "potentials/patchyParticleSTMV.hpp"
#include "potentials/patchyProtein.hpp"
#include "potentials/patchyProteinMarkovSwitch.hpp"
#include "potentials/patchyProteinMAPK.hpp"


namespace msmrd {
    /*
     * pyBinders for the c++ potentials classes (see bindInternal for the parent classes)
     */
    void bindPotentials(py::module &m) {
        pybind11::class_<combinedPairPotential, pairPotential>(m, "combinedPairPotential")
                .def(py::init<>())
                .def("addPotential", &combinedPairPotential::addPotential);

        pybind11::class_<combinedExternalPotential, externalPotential>(m, "combinedExternalPotential")
                .def(py::init<>())
                .def("addPotential", &combinedExternalPotential::addPotential);

        py::class_<dipole, externalPotential>(m, "dipole", "dipole potential (scalefactor, "
                                                           "directionEField)")
                .def(py::init<double &, std::vector<double> &>());

        py::class_<gaussians3D, externalPotential >(m, "gaussians3D", "Gaussians-based potential (nminima, "
                                                                        "maxrad, scalefactor, seed), alternative"
                                                                        "constructors available.")
                .def(py::init<int &, double &, double &, long &>())
                .def(py::init<std::vector<std::vector<double>> &,
                        std::vector<std::vector<double>>&, double &>())
                .def(py::init<std::vector<std::vector<double>> &,
                        std::vector<std::vector<double>>&, std::vector<int> &, double &>());

        py::class_<harmonic, externalPotential >(m, "harmonic", "Harmonic potential (minima, "
                                                                      "konstants, scalefactor), alternative"
                                                                      "constructors available.")
                .def(py::init<std::vector<double> &,std::vector<double>&, double &>())
                .def(py::init<std::vector<double> &,std::vector<double>&, std::vector<int> &, double &>());

        py::class_<bistable, externalPotential >(m, "bistable", "Bistable potential (minimaDist, "
                                                                "konstants, scalefactor), alternative"
                                                                "constructors available.")
                .def(py::init<double &,std::vector<double>&, double &>())
                .def(py::init<double &,std::vector<double>&, std::vector<int> &, double &>());

        py::class_<bistable2, externalPotential >(m, "bistable2", "Bistable potential 2 (parameters, "
                                                                "scalefactor), alternative"
                                                                "bistable potential.")
                .def(py::init<std::vector<double>&, double &>())
                .def(py::init<std::vector<double>&, std::vector<int> &, double &>());

        py::class_<lennardJones, pairPotential>(m, "lennardJones", "Lennard-Jones potential "
                                                           "(epsilon, sigma)")
                .def(py::init<double &, double &>())
                .def(py::init<double &, double &, double &>())
                .def(py::init<double &, double &, std::vector<int> &>())
                .def("setForceCapValue", &lennardJones::setForceCapValue)
                .def("setPotentialCutOff", &lennardJones::setPotentialCutOff)
                .def("getPotentialCutOff", &lennardJones::getPotentialCutOff)
                .def("getForceCapValue", &lennardJones::getForceCapValue);

        py::class_<pairBistable, pairPotential>(m, "pairBistable", "pairBistable potential "
                                                                   "(x0, rad, scalefactor)")
                .def(py::init<double &, double &, double &>())
                .def(py::init<double &, double &, std::vector<int> &, double &>());

        py::class_<WCA, lennardJones>(m, "WCA", "WCA potential (epsilon, sigma)")
                .def(py::init<double &, double &>())
                .def(py::init<double &, double &, std::vector<int> &>());

        py::class_<harmonicRepulsion, pairPotential>(m, "HarmonicRepulsion", "Harmonic repulsion potential "
                                                                   "(k, range)")
                .def(py::init<double &, double &>());

        py::class_<gayBerne, pairPotential>(m, "gayBerne", "Gay-Berne potential "
                                                           "(a, d, eps0, sig0)")
                .def(py::init<double &, double &, double &, double &>());


        py::class_<patchyParticle, pairPotential>(m, "patchyParticle", "Patchy "
                                                     "particle potential (sigma, strength, patchesCoordinates)")
                .def(py::init<double &, double &, std::vector<std::vector<double>> &>());


        py::class_<patchyParticleAngular, patchyParticle>(m, "patchyParticleAngular", "Patchy particle potential with"
                                                             "explicit angular dependence(sigma, strength, "
                                                             "patchesCoordinates)")
                .def(py::init<double &, double &, std::vector<std::vector<double>> &>())
                .def(py::init<double &, double &, double &, std::vector<std::vector<double>> &>());


        py::class_<patchyParticleAngular2, patchyParticleAngular>(m, "patchyParticleAngular2", "Patchy particle "
                                                                     "potential with explicit angular dependence in "
                                                                     "terms of quaternions (sigma, strength, "
                                                                     "patchesCoordinates)")
                .def(py::init<double &, double &, std::vector<std::vector<double>> &>())
                .def(py::init<double &, double &, double &, std::vector<std::vector<double>> &>());

        py::class_<patchyParticleSTMV, pairPotential>(m, "patchyParticleSTMV", "Patchy "
                                                      "particle potential for Satellite Tobacco"
                                                      "Mosaic virus (sigma, strength, angularStrength)")
                .def(py::init<double &, double &>())
                .def("getPartPosition", &patchyParticleSTMV::getPartPosition)
                .def("setParticlesDiameters", &patchyParticleSTMV::setParticlesDiameters)
                .def("setRepulsivePotentialParameters", &patchyParticleSTMV::setRepulsivePotentialParameters)
                .def("setInteractingPatchesPotentialParameters",
                     &patchyParticleSTMV::setInteractingPatchesPotentialParameters);

        py::class_<patchyProtein, pairPotential>(m, "patchyProtein",
                                                 "Patchy protein potential (sigma, strength, patches"
                                                 " coordinates A, patches coordinates B)")
                .def(py::init<double &, double &,
                        std::vector<std::vector<double>> &,
                        std::vector<std::vector<double>> &>())
                .def("setAttractivePotentialParameters", &patchyProtein::setAttractivePotentialParameters)
                .def("setRepulsivePotentialParameters", &patchyProtein::setRepulsivePotentialParameters)
                .def("setInteractingPatchesPotentialParameters",
                        &patchyProtein::setInteractingPatchesPotentialParameters);


        py::class_<patchyProteinMarkovSwitch, patchyProtein>(m, "patchyProteinMarkovSwitch",
                                                 "Patchy protein with Markov Switch potential (sigma, "
                                                 "strength, angular strength, patches coordinates A, patches coordinates B)")
                .def(py::init<double &, double &, double &,
                        std::vector<std::vector<double>> &,
                        std::vector<std::vector<double>> &>());

        py::class_<patchyProteinMAPK, patchyProtein>(m, "patchyProteinMAPK", "patchy protein potenial for MAPK "
                                                        "simulation (sigma, strength, patches coordinates A, "
                                                        "patches coordinates B)")
                .def(py::init<double &, double &,
                        std::vector<std::vector<double>> &,
                        std::vector<std::vector<double>> &>());

    }

}
