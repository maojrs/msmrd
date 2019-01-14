#include "binding.hpp"
#include "trajectories/trajectory.hpp"


namespace msmrd {
    /*
     * pyBinders for the c++ trajectories classes
     */
    void bindTrajectories(py::module &m) {

        py::class_<trajectoryPosition, trajectory>(m, "trajectoryPosition", "position trajectory (#particles or #pairs of particles, "
                                                                "approx size)")
                .def(py::init<int &, int &>())
                .def("sample", &trajectoryPosition::sample)
                .def("sampleRelative", &trajectoryPosition::sampleRelative)
                .def("write2file", &trajectoryPosition::write2file<4>)
                .def("write2H5file", &trajectoryPosition::write2H5file<4>)
                .def("writeChunk2H5file", &trajectoryPosition::writeChunk2H5file<4>)
                .def("emptyBuffer", &trajectoryPosition::emptyBuffer)
                .def_property_readonly("data", &trajectoryPosition::getData);


        py::class_<trajectoryPositionOrientation, trajectory>(m, "trajectoryPositionOrientation", "position and orientation "
                                                                                      "trajectory (#particles or #pairs "
                                                                                      "of particles, approx size)")
                .def(py::init<int &, int &>())
                .def("sample", &trajectoryPositionOrientation::sample)
                .def("sampleRelative", &trajectoryPositionOrientation::sampleRelative)
                .def("write2file", &trajectoryPositionOrientation::write2file<8>)
                .def("write2H5file", &trajectoryPositionOrientation::write2H5file<8>)
                .def("writeChunk2H5file", &trajectoryPositionOrientation::writeChunk2H5file<8>)
                .def("emptyBuffer", &trajectoryPositionOrientation::emptyBuffer)
                .def_property_readonly("data", &trajectoryPositionOrientation::getData);
    }

}