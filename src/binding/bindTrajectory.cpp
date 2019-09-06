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
                .def("write2H5file", &trajectoryPosition::write2H5file<double, 4>)
                .def("writeChunk2H5file", &trajectoryPosition::writeChunk2H5file<double, 4>);


        py::class_<trajectoryPositionOrientation, trajectory>(m, "trajectoryPositionOrientation", "position and "
                                                                                                  "orientation trajectory "
                                                                                                  "(#particles or #pairs "
                                                                                                  "of particles, approx size)")
                .def(py::init<int &, int &>())
                .def("write2H5file", &trajectoryPositionOrientation::write2H5file<double, 8>)
                .def("writeChunk2H5file", &trajectoryPositionOrientation::writeChunk2H5file<double, 8>);



        /* Not defined as child class since parent class is a virtual template, so need to add all functions
         * from parent classes manually (all the way to the original trajectory.hpp parent)*/
        py::class_<patchyDimer>(m, "patchyDimer", "discrete trajectory of patchy dimer"
                                                  "example(#particles or #pairs of particles, approx size)")
                .def(py::init<int &, int &>())
                .def(py::init<int &, int &, double &, double &>())
                .def("setBoundary", &patchyDimer::setBoundary)
                .def("sample", &patchyDimer::sample)
                .def("sampleRelative", &patchyDimer::sampleRelative)
                .def("write2file", &patchyDimer::write2file<double>)
                .def("emptyBuffer", &patchyDimer::emptyBuffer)
                .def("sampleDiscreteTrajectory", &patchyDimer::sampleDiscreteTrajectory)
                .def("getState", &patchyDimer::getState)
                .def("discretizeTrajectory", &patchyDimer::discretizeTrajectory)
                .def("discretizeTrajectoryH5", &patchyDimer::discretizeTrajectoryH5)
                .def("write2H5file", &patchyDimer::write2H5file<double, 8>)
                .def("writeChunk2H5file", &patchyDimer::writeChunk2H5file<double, 8>);

        // Not defined as child class since prent class is a virtual template
        py::class_<patchyDimer2>(m, "patchyDimer2", "discrete trajectory of patchy dimer 2"
                                                  "example(#particles or #pairs of particles, approx size)")
                .def(py::init<int &, int &>())
                .def(py::init<int &, int &, double &, double &>())
                .def("setBoundary", &patchyDimer2::setBoundary)
                .def("sample", &patchyDimer2::sample)
                .def("sampleRelative", &patchyDimer2::sampleRelative)
                .def("write2file", &patchyDimer2::write2file<double>)
                .def("emptyBuffer", &patchyDimer2::emptyBuffer)
                .def("sampleDiscreteTrajectory", &patchyDimer2::sampleDiscreteTrajectory)
                .def("getState", &patchyDimer2::getState)
                .def("discretizeTrajectory", &patchyDimer2::discretizeTrajectory)
                .def("discretizeTrajectoryH5", &patchyDimer2::discretizeTrajectoryH5)
                .def("write2H5file", &patchyDimer2::write2H5file<double, 8>)
                .def("writeChunk2H5file", &patchyDimer2::writeChunk2H5file<double, 8>);



    }
}