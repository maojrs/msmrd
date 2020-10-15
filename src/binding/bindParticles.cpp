//
// Created by maojrs on 9/7/18.
//
#include "binding.hpp"
#include "particle.hpp"
#include "particleCompound.hpp"


namespace msmrd {
    /*
     * pyBinders for the c++ particles classes
     */
    void bindParticles(py::module &m) {
        py::class_<particle>(m, "particle", "particle (D, Drot, position, orientation (quaternion)) or"
                                            "(type, state, D, Drot, position, orientation (quaternion)) "
                                            "for MArkovian switch behavior")
                .def(py::init<double &, double &, std::vector<double> &, std::vector<double> &>())
                .def(py::init<int &, int &, double &, double &, std::vector<double> &, std::vector<double> &>())
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
                .def_property_readonly("nextPosition", [](const particle &part) {
                    return vec2numpy(3, part.nextPosition);
                }, "particle position")
                .def_property_readonly("nextOrientation", [](const particle &part) {
                    return vec2numpy(4, part.nextOrientation);
                })
                .def_property_readonly("state", &particle::getState)
                .def_property_readonly("lagtime", &particle::getLagtime)
                .def_property_readonly("isMSMactive", &particle::isMSMactive)
                .def_property_readonly("boundTo", &particle::getBoundTo)
                .def_property_readonly("boundState", &particle::getBoundState)
                .def_property_readonly("boundList", &particle::getBoundList)
                .def_property_readonly("boundStates", &particle::getBoundStates)
                .def_property_readonly("compoundIndex", &particle::getCompoundIndex)


                .def("activate", &particle::activate)
                .def("deactivate", &particle::deactivate)
                .def("setD", &particle::setD)
                .def("setDrot", &particle::setDrot)
                .def("setType", &particle::setType)
                .def("setPosition", &particle::setPositionPyBind)
                .def("setOrientVector", &particle::setOrientVectorPyBind)
                .def("setOrientation", &particle::setOrientationPyBind)
                .def("setBoundTo", &particle::setBoundTo)
                .def("setBoundState", &particle::setBoundState)
                .def("activateMSM", &particle::activateMSM)
                .def("deactivateMSM", &particle::deactivateMSM)
                .def("setState", &particle::setState)
                .def("setLagtime", &particle::setLagtime)
                .def("setMSMoff", &particle::setMSMoff)
                .def("setMSMon", &particle::setMSMon)
                .def("setActivePatchList", &particle::setActivePatchList);


        py::class_<particleCompound>(m, "particleCompound", "particle complex class that keeps track of"
                                                          "all bound complexes in multiparticle MSM/RD.")
                .def(py::init<>())
                .def(py::init<std::vector<double> &>())
                .def(py::init<std::map<std::tuple<int,int>, int> &>())
                .def(py::init<std::vector<double> &, std::map<std::tuple<int,int>, int> &>());
    }

}


