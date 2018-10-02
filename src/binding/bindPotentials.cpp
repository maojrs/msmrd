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
     * pyBinders for the c++ potentials classes
     */
    void bindPotentials(py::module &m) {
        /* Bind external potential parent classes as instantiated templates
         * so inheritance works correctly. */
        pybind11::class_<externalPotential<> >(m, "externalPotential");
        pybind11::class_<externalPotential<int> >(m, "externalMixPotential");
        pybind11::class_<externalPotential<vec3<double>>>(m, "externalRodPotential");
        pybind11::class_<externalPotential<vec3<double>, int>>(m, "externalRodMixPotential");
        pybind11::class_<externalPotential<quaternion<double>>>(m, "externalRigidBodyPotential");
        pybind11::class_<externalPotential<quaternion<double>, int>>(m, "externalRigidBodyMixPotential");

        /* Bind pair potential parent classes as instantiated templates
         * so inheritance works correctly. */
        pybind11::class_<pairPotential<> >(m, "pairPotential");
        pybind11::class_<pairPotential<int, int> >(m, "pairMixPotential");
        pybind11::class_<pairPotential<vec3<double>, vec3<double>>>(m, "pairRodPotential");
        pybind11::class_<pairPotential<vec3<double>, vec3<double>, int, int>>(m, "pairRodMixPotential");
        pybind11::class_<pairPotential<quaternion<double>, quaternion<double>>>(m, "pairRigidBodyPotential");
        pybind11::class_<pairPotential<quaternion<double>, quaternion<double>, int, int>>
                (m, "pairRigidBodyMixPotential");


        py::class_<gaussians3D, externalPotential<> >(m, "gaussians3D")
                .def(py::init<int &, double &, double &, long &>())
                .def("evaluate", (double (gaussians3D::*)
                        (std::vector<double>)) &gaussians3D::evaluatePyBind)
                .def("forceTorque", (std::vector<std::vector<double>> (gaussians3D::*)
                        (std::vector<double>))&gaussians3D::forceTorquePyBind);

        py::class_<dipole, externalPotential<vec3<double>>>(m, "dipole")
                .def(py::init<double &, std::vector<double> &>())
                .def("evaluate", (double (dipole::*)
                        (std::vector<double>, std::vector<double>)) &dipole::evaluatePyBind)
                .def("forceTorque", (std::vector<std::vector<double>> (dipole::*)
                        (std::vector<double>, std::vector<double>)) &dipole::forceTorquePyBind);

        py::class_<gayBerne, pairPotential<vec3<double>, vec3<double>>>(m, "gayBerne")
                .def(py::init<double &, double &, double &, double &>())
                .def("evaluate", (double (gayBerne::*)
                        (std::vector<double>, std::vector<double>,
                         std::vector<double>, std::vector<double>))
                        &gayBerne::evaluatePyBind)
                .def("forceTorque", (std::vector<std::vector<double>> (gayBerne::*)
                        (std::vector<double>, std::vector<double>,
                         std::vector<double>, std::vector<double>))
                        &gayBerne::forceTorquePyBind);

        py::class_<patchyParticle, pairPotential<quaternion<double>, quaternion<double>>>(m, "patchyParticle")
                .def(py::init<double &, double &, std::vector<std::vector<double>> &>())
                .def("evaluate", (double (patchyParticle::*)
                        (std::vector<double>, std::vector<double>,
                         std::vector<double>, std::vector<double>))
                        &patchyParticle::evaluatePyBind)
                .def("forceTorque", (std::vector<std::vector<double>> (patchyParticle::*)
                        (std::vector<double>, std::vector<double>,
                         std::vector<double>, std::vector<double>))
                        &patchyParticle::forceTorquePyBind);

//        py::class_<patchyProtein, pairPotential<quaternion<double>, quaternion<double>, int, int>>(m, "patchyProtein")
//                .def(py::init<double &, double &, std::vector<std::vector<double>> &>())
//                .def("evaluate", (double (patchyParticle::*)
//                        (std::vector<double>, std::vector<double>,
//                         std::vector<double>, std::vector<double>))
//                        &patchyParticle::evaluatePyBind)
//                .def("forceTorque", (std::vector<std::vector<double>> (patchyParticle::*)
//                        (std::vector<double>, std::vector<double>,
//                         std::vector<double>, std::vector<double>))
//                        &patchyParticle::forceTorquePyBind);
    }

}
