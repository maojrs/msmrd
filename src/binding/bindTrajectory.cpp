#include "binding.hpp"
#include "trajectories/trajectoryPosition.hpp"
#include "trajectories/trajectoryPositionOrientation.hpp"


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
                .def("write2file", &trajectoryPosition::write2file<double>)
                .def("write2H5file", &trajectoryPosition::write2H5file<double,4>)
                .def("writeChunk2H5file", &trajectoryPosition::writeChunk2H5file<double,4>)
                .def("emptyBuffer", &trajectoryPosition::emptyBuffer)
                .def_property_readonly("data", &trajectoryPosition::getTrajectoryData);


        py::class_<trajectoryPositionOrientation, trajectory>(m, "trajectoryPositionOrientation", "position and orientation "
                                                                                      "trajectory (#particles or #pairs "
                                                                                      "of particles, approx size)")
                .def(py::init<int &, int &>())
                .def("sample", &trajectoryPositionOrientation::sample)
                .def("sampleRelative", &trajectoryPositionOrientation::sampleRelative)
                .def("write2file", &trajectoryPositionOrientation::write2file<double>)
                .def("write2H5file", &trajectoryPositionOrientation::write2H5file<double,8>)
                .def("writeChunk2H5file", &trajectoryPositionOrientation::writeChunk2H5file<double,8>)
                .def("emptyBuffer", &trajectoryPositionOrientation::emptyBuffer)
                .def_property_readonly("data", &trajectoryPositionOrientation::getTrajectoryData);
    }

}