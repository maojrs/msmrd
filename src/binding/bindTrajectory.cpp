#include "binding.hpp"
#include "trajectories/trajectory.hpp"


namespace msmrd {
    /*
     * pyBinders for the c++ trajectories classes
     */
    void bindTrajectories(py::module &m) {

        py::class_<trajectoryPosition>(m, "trajectoryPositionOrientation")
                .def(py::init<int &, int &>())
                .def("sample", &trajectoryPosition::sample);

        py::class_<trajectoryPositionOrientation>(m, "trajectoryPositionOrientation")
                .def(py::init<int &, int &>())
                .def("sample", &trajectoryPositionOrientation::sample);

        py::class_<twoParticleRelativeTrajectory>(m, "twoParticleRelativeTrajectory")
                .def(py::init<int &>())
                .def("sample", &twoParticleRelativeTrajectory::sample);
    }

}