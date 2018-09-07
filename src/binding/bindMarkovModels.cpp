//
// Created by maojrs on 9/7/18.
//
#include "binding.hpp"
#include "msm.hpp"

namespace msmrd {
    using msm = msmrd::discreteTimeMarkovStateModel;
    using ctmsm = msmrd:: continuousTimeMarkovStateModel;
    /*
     * pyBinders for the c++ Markov state models (MSMs) classes
     */
    void bindMarkovModels(py::module &m) {
        py::class_<msm>(m, "msm")
                .def(py::init<int &, std::vector<std::vector<double>> &, double &, long &>())
                .def_property_readonly("ID", &msm::getID)
                .def_property_readonly("nstates", &msm::getNstates)
                .def_property_readonly("lagtime", &msm::getLagtime)
                .def_property_readonly("tmatrix", &msm::getTmatrix)
                .def_property_readonly("D", [](const msm &currentmsm) {
                    return vec2numpy(currentmsm.nstates, currentmsm.Dlist);
                })
                .def_property_readonly("Drot", [](const msm &currentmsm) {
                    return vec2numpy(currentmsm.nstates, currentmsm.Drotlist);
                })
                .def("setD", &msm::setD)
                .def("setDrot", &msm::setDrot)
                .def("propagate", &msm::propagate);

        py::class_<ctmsm>(m, "ctmsm")
                .def(py::init<int &, std::vector<std::vector<double>> &, long &>())
                .def_property_readonly("ID", &ctmsm::getID)
                .def_property_readonly("nstates", &ctmsm::getNstates)
                .def_property_readonly("lagtime", &ctmsm::getLagtime)
                .def_property_readonly("tmatrix", &ctmsm::getTmatrix)
                .def_property_readonly("D", [](const ctmsm &currentmsm) {
                    return vec2numpy(currentmsm.nstates, currentmsm.Dlist);
                })
                .def_property_readonly("Drot", [](const msm &currentmsm) {
                    return vec2numpy(currentmsm.nstates, currentmsm.Drotlist);
                })
                .def("setD", &ctmsm::setD)
                .def("setDrot", &ctmsm::setDrot)
                .def("propagate", &ctmsm::propagate);
    }

}
