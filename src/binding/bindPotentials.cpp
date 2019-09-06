#include "binding.hpp"
#include "potentials/potentials.hpp"
#include "potentials/dipole.hpp"
#include "potentials/gaussians3D.hpp"
#include "potentials/harmonicRepulsion.hpp"
#include "potentials/gayBerne.hpp"
#include "potentials/patchyParticle.hpp"
#include "potentials/patchyParticleAngular.hpp"
#include "potentials/patchyParticleAngular2.hpp"
#include "potentials/patchyProtein.hpp"
#include "potentials/patchyProteinMarkovSwitch.hpp"


namespace msmrd {
    /*
     * pyBinders for the c++ potentials classes (see bindInternal for the parent classes)
     */
    void bindPotentials(py::module &m) {
        py::class_<gaussians3D, externalPotential >(m, "gaussians3D", "Gaussians-based potential (nminima, "
                                                                        "maxrad, scalefactor, seed)")
                .def(py::init<int &, double &, double &, long &>());


        py::class_<dipole, externalPotential>(m, "dipole", "dipole potential (scalefactor, "
                                                                         "directionEField)")
                .def(py::init<double &, std::vector<double> &>());


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


        py::class_<patchyProtein, pairPotential>(m, "patchyProtein",
                                                 "Patchy protein potential (sigma, strength, patches"
                                                 " coordinates A, patches coordinates B)")
                .def(py::init<double &, double &,
                        std::vector<std::vector<double>> &,
                        std::vector<std::vector<double>> &>());


        py::class_<patchyProteinMarkovSwitch, patchyProtein>(m, "patchyProteinMarkovSwitch",
                                                 "Patchy protein with Markov Switch potential (sigma, "
                                                 "strength, patches coordinates A, patches coordinates B)")
                .def(py::init<double &, double &,
                        std::vector<std::vector<double>> &,
                        std::vector<std::vector<double>> &>());
    }

}
