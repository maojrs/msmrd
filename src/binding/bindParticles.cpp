//
// Created by maojrs on 9/7/18.
//
#include "binding.hpp"
#include "particle.hpp"

namespace msmrd {
    /*
     * pyBinders for the c++ particles classes
     */
    void bindParticles(py::module &m) {
        py::class_<particle>(m, "particle", "particle (D, Drot, position, orientation (quaternion))")
                .def(py::init<double &, double &, std::vector<double> &, std::vector<double> &>())
                .def_property_readonly("ID", &particle::getID)
                .def_property_readonly("D", &particle::getD)
                .def_property_readonly("Drot", &particle::getDrot)
                .def_property_readonly("type", &particle::getType)
                .def_property_readonly("isActive", &particle::isActive)
                .def_property_readonly("position", [](const particle &part) {
                    return vec2numpy(3, part.position);
                }, "particle position")
                .def_property_readonly("orientvector", [](const particle &part) {
                    return vec2numpy(3, part.orientvector);
                })
                .def_property_readonly("orientation", [](const particle &part) {
                    return vec2numpy(4, part.orientation);
                })
                .def("activate", &particle::activate)
                .def("deactivate", &particle::deactivate)
                .def("setD", &particle::setD)
                .def("setDrot", &particle::setDrot)
                .def("setType", &particleMS::setType)
                .def("setPosition", &particle::setPositionPyBind)
                .def("setOrientVector", &particle::setOrientVectorPyBind)
                .def("setOrientation", &particle::setOrientationPyBind)
                .def("setSymmetryQuaternions", &particle::setSymmetryQuaternions);


        py::class_<particleMS, particle>(m, "particleMS", "particle with Markovian switch (type, state, "
                                                          "D, Drot, position, orientation (quaternion))")
                .def(py::init<int &, int &, double &, double &, std::vector<double> &, std::vector<double> &>())
                .def_property_readonly("state", &particleMS::getState)
                .def_property_readonly("lagtime", &particleMS::getLagtime)
                .def_property_readonly("isMSMactive", &particleMS::isMSMactive)
                .def("activateMSM", &particleMS::activateMSM)
                .def("deactivateMSM", &particleMS::deactivateMSM)
                .def("setState", &particleMS::setState)
                .def("setLagtime", &particleMS::setLagtime)
                .def("setMSMoff", &particleMS::setMSMoff)
                .def("setMSMon", &particleMS::setMSMon);
    }

}


