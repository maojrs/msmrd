#include "binding.hpp"
#include "potentials/potentials.hpp"
#include "potentials/dipole.hpp"
#include "potentials/gaussians3D.hpp"
#include "potentials/harmonicRepulsion.hpp"
#include "potentials/gayBerne.hpp"
#include "potentials/patchyParticle.hpp"
#include "potentials/patchyProtein.hpp"


namespace msmrd {
    /*
     * pyBinders for the c++ potentials classes (see bindInternal for the parent classes)
     */
    void bindPotentials(py::module &m) {
        py::class_<gaussians3D, externalPotential >(m, "gaussians3D", "Gaussians-based potential (nminima, "
                                                                        "maxrad, scalefactor, seed)")
                .def(py::init<int &, double &, double &, long &>())
                .def("evaluate", &gaussians3D::evaluate)
                .def("forceTorque", &gaussians3D::forceTorque);

        py::class_<dipole, externalPotential>(m, "dipole", "dipole potential (scalefactor, "
                                                                         "directionEField)")
                .def(py::init<double &, std::vector<double> &>())
                .def("evaluate", &dipole::evaluate)
                .def("forceTorque", &dipole::forceTorque);

        py::class_<gayBerne, pairPotential>(m, "gayBerne", "Gay-Berne potential "
                                                           "(a, d, eps0, sig0)")
                .def(py::init<double &, double &, double &, double &>())
                .def("evaluate", &gayBerne::evaluate)
                .def("forceTorque", &gayBerne::forceTorque);

        py::class_<patchyParticle, pairPotential>(m, "patchyParticle", "Patchy "
                                                     "particle potential (sigma, strength, patchesCoordinates)")
                .def(py::init<double &, double &, std::vector<std::vector<double>> &>())
                .def("evaluate", &patchyParticle::evaluate)
                .def("forceTorque", &patchyParticle::forceTorque);

        py::class_<patchyProtein, pairPotential>(m, "patchyProtein",
                                                 "Patchy protein potential (sigma, strength, patches"
                                                 " coordinates A, patches coordinates B)")
                .def(py::init<double &, double &,
                        std::vector<std::vector<double>> &,
                        std::vector<std::vector<double>> &>())
                .def("evaluate", &patchyProtein::evaluate)
                .def("forceTorque", &patchyProtein::forceTorque);
    }

}
