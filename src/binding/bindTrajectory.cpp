#include "binding.hpp"
#include "trajectories/trajectory.hpp"


namespace msmrd {
    /*
     * pyBinders for the c++ trajectories classes
     */
    void bindTrajectories(py::module &m) {

        py::class_<trajectoryPosition>(m, "trajectoryPosition", "position trajectory (#particles or #pairs of particles, "
                                                                "approx size)")
                .def(py::init<int &, int &>())
                .def("sample", &trajectoryPosition::sample)
                .def("sampleRelative", &trajectoryPosition::sampleRelative)
                .def("write2file", &trajectoryPosition::write2file)
                .def_property_readonly("data", &trajectoryPosition::getData);


        py::class_<trajectoryPositionOrientation>(m, "trajectoryPositionOrientation", "position and orientation "
                                                                                      "trajectory (#particles or #pairs "
                                                                                      "of particles, approx size)")
                .def(py::init<int &, int &>())
                .def("sample", &trajectoryPositionOrientation::sample)
                .def("sampleRelative", &trajectoryPositionOrientation::sampleRelative)
                .def_property_readonly("data", &trajectoryPositionOrientation::getData);
    }

}