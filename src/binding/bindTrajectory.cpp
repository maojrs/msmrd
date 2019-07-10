#include "binding.hpp"
#include "trajectories/trajectoryPosition.hpp"
#include "trajectories/trajectoryPositionOrientation.hpp"
#include "trajectories/discrete/patchyDimer.hpp"



namespace msmrd {
    /*
     * pyBinders for the c++ trajectories classes
     */
    void bindTrajectories(py::module &m) {

        py::class_<trajectoryPosition, trajectory>(m, "trajectoryPosition", "position trajectory (#particles or "
                                                                            "#pairs of particles, approx size)")
                .def(py::init<int &, int &>())
                .def("write2H5file", &trajectoryPosition::write2H5file<double,4>)
                .def("writeChunk2H5file", &trajectoryPosition::writeChunk2H5file<double,4>);



        py::class_<trajectoryPositionOrientation, trajectory>(m, "trajectoryPositionOrientation", "position and "
                                                                                      "orientation trajectory "
                                                                                      "(#particles or #pairs "
                                                                                      "of particles, approx size)")
                .def(py::init<int &, int &>())
                .def("write2H5file", &trajectoryPositionOrientation::write2H5file<double,8>)
                .def("writeChunk2H5file", &trajectoryPositionOrientation::writeChunk2H5file<double,8>);


        py::class_<patchyDimer, trajectoryPositionOrientation>(m, "patchyDimer", "discrete trajectory of patchy dimer"
                                                                                 "exampole(#particles or #pairs "
                                                                                 "of particles, approx size)")
                .def(py::init<int &, int &>())
                .def("sampleDiscreteTrajectory", &patchyDimer::sampleDiscreteTrajectory)
                .def("getState", &patchyDimer::getState)
                .def("discretizeTrajectory", &patchyDimer::discretizeTrajectory);


    }

}