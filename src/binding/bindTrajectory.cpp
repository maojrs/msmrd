#include "binding.hpp"
#include "trajectories/trajectory.hpp"


namespace msmrd {
    /*
     * pyBinders for the c++ integrators classes
     */
    void bindSimulation(py::module &m) {
        py::class_<trajectoryPositionOrientation>(m, "trajectoryPositionOrientation")
                .def(py::init<int &, int &>());
                //.def("get_sampler", );

    }

}