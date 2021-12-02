//
// Created by maojrs on 10/8/18.
//
#include "binding.hpp"
#include "boundaries/boundary.hpp"
#include "integrators/integrator.hpp"
#include "integrators/integratorLangevin.hpp"
#include "markovModels/markovModel.hpp"
#include "potentials/potentials.hpp"
#include "trajectories/trajectory.hpp"



namespace msmrd {
    /*
     * pyBinders for the c++ parent classes (mostly virtual) that should remain hidden under _internal in python
     */
    void bindInternal(py::module &m) {

        /* Bind boundaries parent class*/
        pybind11::class_<boundary>(m, "boundary")
                .def_property_readonly("boundaryType", &boundary::getBoundaryType);

        /* Bind integrators parent class*/
        pybind11::class_<integrator>(m, "integrator")
                .def_property_readonly("clock", &integrator::getClock)
                .def("setClock", &integrator::setClock)
                .def("resetClock", &integrator::resetClock)
                .def("setKbT", &integrator::setKbT)
                .def("setBoundary", &integrator::setBoundary)
                .def("setExternalPotential", &integrator::setExternalPotential)
                .def("setPairPotential", &integrator::setPairPotential);

        /* Bind Langevin integrators parent class*/
        pybind11::class_<integratorLangevin, integrator>(m, "integratorLangevin")
                .def_property_readonly("clock", &integratorLangevin::getClock)
                .def("setClock", &integratorLangevin::setClock)
                .def("resetClock", &integratorLangevin::resetClock)
                .def("setKbT", &integratorLangevin::setKbT)
                .def("setBoundary", &integratorLangevin::setBoundary)
                .def("setExternalPotential", &integratorLangevin::setExternalPotential)
                .def("setPairPotential", &integratorLangevin::setPairPotential);

        /* Bind Markov models parent class*/
        pybind11::class_<markovModel>(m, "markovModel")
                .def_property_readonly("ID", &markovModel::getID)
                .def_property_readonly("nstates", &markovModel::getNstates)
                .def_property_readonly("lagtime", &markovModel::getLagtime)
                .def_property_readonly("tmatrix", &markovModel::getTmatrix)
                .def_property_readonly("D", [](const markovModel &currentmsm) {
                    return vec2numpy(currentmsm.nstates, currentmsm.Dlist);
                })
                .def_property_readonly("Drot", [](const markovModel &currentmsm) {
                    return vec2numpy(currentmsm.nstates, currentmsm.Drotlist);
                })
                .def("setD", &markovModel::setD)
                .def("setDrot", &markovModel::setDrot)
                .def("getTransitionMatrix", &markovModel::getTmatrix);

        /* Bind trajectories parent class*/
        pybind11::class_<trajectory>(m, "trajectory")
                .def_property_readonly("data", &trajectory::getTrajectoryData)
                .def("setBoundary", &trajectory::setBoundary)
                .def("sample", &trajectory::sample)
                .def("sampleRelative", &trajectory::sampleRelative)
                .def("write2file", &trajectory::write2file<double>)
                .def("emptyBuffer", &trajectory::emptyBuffer);

        /* Bind external potential parent class  */
        pybind11::class_<externalPotential>(m, "externalPotential")
                .def("evaluate", &externalPotential::evaluate)
                .def("forceTorque", &externalPotential::forceTorquePyBind);

        /* Bind pair potential parent class */
        pybind11::class_<pairPotential>(m, "pairPotential")
                .def("evaluate", &pairPotential::evaluate)
                .def("forceTorque", &pairPotential::forceTorquePyBind);
    }

}
