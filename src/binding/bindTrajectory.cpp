#include "binding.hpp"
#include "trajectories/trajectoryEnergyTemperature.hpp"
#include "trajectories/trajectoryPosition.hpp"
#include "trajectories/trajectoryPositionOrientation.hpp"
#include "trajectories/trajectoryPositionVelocity.hpp"
#include "trajectories/discrete/patchyDimerTrajectory.hpp"
#include "trajectories/discrete/patchyProteinTrajectory.hpp"
#include "trajectories/discrete/MAPKtrajectory.hpp"


namespace msmrd {
    /*
     * pyBinders for the c++ trajectories classes
     */
    void bindTrajectories(py::module &m) {

        py::class_<trajectoryEnergyTemperature, trajectory>(m, "trajectoryEnergyTemperature", "energy temperature trajectory"
                                                                                              "(approx size)")
                .def(py::init<int &>())
                .def("write2H5file", &trajectoryEnergyTemperature::write2H5file<double, 3>)
                .def("writeChunk2H5file", &trajectoryEnergyTemperature::writeChunk2H5file<double, 3>);

        py::class_<trajectoryPosition, trajectory>(m, "trajectoryPosition", "position trajectory (#particles or "
                                                                            "#pairs of particles, approx size)")
                .def(py::init<int &, int &>())
                .def("write2H5file", &trajectoryPosition::write2H5file<double, 4>)
                .def("writeChunk2H5file", &trajectoryPosition::writeChunk2H5file<double, 4>);

        py::class_<trajectoryPositionType, trajectoryPosition>(m, "trajectoryPositionType",
                "position type trajectory (#particles or #pairs of particles, approx size)")
                .def(py::init<int &, int &>())
                .def("write2H5file", &trajectoryPositionType::write2H5file<double, 5>)
                .def("writeChunk2H5file", &trajectoryPositionType::writeChunk2H5file<double, 5>);

        py::class_<trajectoryPositionDistinguished, trajectoryPositionType>(m, "trajectoryPositionDistinguished",
                                                               "position trajectory for distinguished particle "
                                                               "(#particles or #pairs of particles, approx size, "
                                                               "distinguihed types)")
                .def(py::init<int &, int &, std::vector<int> &>())
                .def("write2H5file", &trajectoryPositionDistinguished::write2H5file<double, 5>)
                .def("writeChunk2H5file", &trajectoryPositionDistinguished::writeChunk2H5file<double, 5>);

        py::class_<trajectoryPositionOrientation, trajectory>(m, "trajectoryPositionOrientation", "position and "
                                                                                                  "orientation trajectory "
                                                                                                  "(#particles or #pairs "
                                                                                                  "of particles, approx size)")
                .def(py::init<int &, int &>())
                .def("write2H5file", &trajectoryPositionOrientation::write2H5file<double, 8>)
                .def("writeChunk2H5file", &trajectoryPositionOrientation::writeChunk2H5file<double, 8>);

        py::class_<trajectoryPositionOrientationType, trajectoryPositionOrientation>(m,
                "trajectoryPositionOrientationType",
                "position and orientation trajectory (#particles or #pairs of particles, approx size)")
                .def(py::init<int &, int &>())
                .def("write2H5file", &trajectoryPositionOrientationType::write2H5file<double, 9>)
                .def("writeChunk2H5file", &trajectoryPositionOrientationType::writeChunk2H5file<double, 9>);

        py::class_<trajectoryPositionOrientationState, trajectoryPositionOrientation>(m,
                                                      "trajectoryPositionOrientationState", "position, orientation"
                                                      "and state trajectory (#particles or #pairs of particles, "
                                                      "approx size)")
                .def(py::init<int &, int &>())
                .def("write2H5file", &trajectoryPositionOrientationState::write2H5file<double, 9>)
                .def("writeChunk2H5file", &trajectoryPositionOrientationState::writeChunk2H5file<double, 9>);

        py::class_<trajectoryPositionOrientationStateType, trajectoryPositionOrientationState>(m,
                "trajectoryPositionOrientationStateType", "position, orientation, state and type trajectory "
                                                      "(#particles or #pairs of particles, approx size)")
                .def(py::init<int &, int &>())
                .def("write2H5file", &trajectoryPositionOrientationStateType::write2H5file<double, 10>)
                .def("writeChunk2H5file", &trajectoryPositionOrientationStateType::writeChunk2H5file<double, 10>);

        py::class_<trajectoryPositionVelocity, trajectory>(m, "trajectoryPositionVelocity", "position and "
                                                                                            "velocity trajectory "
                                                                                            "(#particles or #pairs "
                                                                                            "of particles, approx size)")
                .def(py::init<int &, int &>())
                .def("write2H5file", &trajectoryPositionVelocity::write2H5file<double, 7>)
                .def("writeChunk2H5file", &trajectoryPositionVelocity::writeChunk2H5file<double, 7>);

        py::class_<trajectoryPositionVelocityType, trajectoryPositionVelocity>(m, "trajectoryPositionVelocityType", "position, "
                                                                                            "velocity and type trajectory "
                                                                                            "(#particles or #pairs "
                                                                                            "of particles, approx size)")
                .def(py::init<int &, int &>())
                .def("write2H5file", &trajectoryPositionVelocityType::write2H5file<double, 8>)
                .def("writeChunk2H5file", &trajectoryPositionVelocityType::writeChunk2H5file<double, 8>);

        py::class_<trajectoryPositionVelocityDistinguished, trajectoryPositionVelocityType>(m, "trajectoryPositionVelocityDistinguished", "position, "
                                                                                               "and velocity trajectory for distinguished"
                                                                                               "particles, approx size, distinguished types)")
                .def(py::init<int &, int &, std::vector<int> &>())
                .def("write2H5file", &trajectoryPositionVelocityDistinguished::write2H5file<double, 8>)
                .def("writeChunk2H5file", &trajectoryPositionVelocityDistinguished::writeChunk2H5file<double, 8>);


        /* Not defined as child class since parent class is a template with virtual fucntions, so need to
         * add all functions from parent classes manually (all the way to the original trajectory.hpp parent).
         * Also note sampleDiscreteState and getState are the same function. */
        py::class_<patchyDimerTrajectory>(m, "patchyDimer", "discrete trajectory of patchy dimer"
                                                  "example(#particles or #pairs of particles, approx size)")
                .def(py::init<int &, int &>())
                .def(py::init<int &, int &, double &, double &>())
                .def("setBoundary", &patchyDimerTrajectory::setBoundary)
                .def("setTolerances", &patchyDimerTrajectory::setTolerances)
                .def("sample", &patchyDimerTrajectory::sample)
                .def("sampleRelative", &patchyDimerTrajectory::sampleRelative)
                .def("write2file", &patchyDimerTrajectory::write2file<double>)
                .def("emptyBuffer", &patchyDimerTrajectory::emptyBuffer)
                .def("sampleDiscreteTrajectory", &patchyDimerTrajectory::sampleDiscreteTrajectory)
                .def("sampleDiscreteState", &patchyDimerTrajectory::sampleDiscreteState)
                .def("getState", &patchyDimerTrajectory::sampleDiscreteState)
                .def("discretizeTrajectory", &patchyDimerTrajectory::discretizeTrajectory)
                .def("discretizeTrajectoryH5", &patchyDimerTrajectory::discretizeTrajectoryH5)
                .def("write2H5file", &patchyDimerTrajectory::write2H5file<double, 8>)
                .def("writeChunk2H5file", &patchyDimerTrajectory::writeChunk2H5file<double, 8>);

        /* Not defined as child class since prent class is a virtual template. Also note sampleDiscreteState and
         * getState are the same function. */
        py::class_<patchyDimerTrajectory2>(m, "patchyDimer2", "discrete trajectory of patchy dimer 2"
                                                  "example(#particles or #pairs of particles, approx size)")
                .def(py::init<int &, int &>())
                .def(py::init<int &, int &, double &, double &>())
                .def("setBoundary", &patchyDimerTrajectory2::setBoundary)
                .def("setTolerances", &patchyDimerTrajectory2::setTolerances)
                .def("sample", &patchyDimerTrajectory2::sample)
                .def("sampleRelative", &patchyDimerTrajectory2::sampleRelative)
                .def("write2file", &patchyDimerTrajectory2::write2file<double>)
                .def("emptyBuffer", &patchyDimerTrajectory2::emptyBuffer)
                .def("sampleDiscreteTrajectory", &patchyDimerTrajectory2::sampleDiscreteTrajectory)
                .def("sampleDiscreteState", &patchyDimerTrajectory2::sampleDiscreteState)
                .def("getState", &patchyDimerTrajectory2::sampleDiscreteState)
                .def("discretizeTrajectory", &patchyDimerTrajectory2::discretizeTrajectory)
                .def("discretizeTrajectoryH5", &patchyDimerTrajectory2::discretizeTrajectoryH5)
                .def("write2H5file", &patchyDimerTrajectory2::write2H5file<double, 8>)
                .def("writeChunk2H5file", &patchyDimerTrajectory2::writeChunk2H5file<double, 8>);


        /* Not defined as child class since parent class is a virtual template, so need to add all functions
         * from parent classes manually (all the way to the original trajectory.hpp parent).
         * Also note sampleDiscreteState and getState are the same function. */
        py::class_<patchyProteinTrajectory>(m, "patchyProtein", "discrete trajectory of patchy protein"
                                                  "example(#particles or #pairs of particles, approx size)")
                .def(py::init<int &, int &>())
                .def(py::init<int &, int &, double &, double &>())
                .def("setBoundary", &patchyProteinTrajectory::setBoundary)
                .def("setTolerances", &patchyProteinTrajectory::setTolerances)
                .def("sample", &patchyProteinTrajectory::sample)
                .def("sampleRelative", &patchyProteinTrajectory::sampleRelative)
                .def("write2file", &patchyProteinTrajectory::write2file<double>)
                .def("emptyBuffer", &patchyProteinTrajectory::emptyBuffer)
                .def("sampleDiscreteTrajectory", &patchyProteinTrajectory::sampleDiscreteTrajectory)
                .def("sampleDiscreteState", &patchyProteinTrajectory::sampleDiscreteState)
                .def("getState", &patchyProteinTrajectory::sampleDiscreteState)
                .def("discretizeTrajectory", &patchyProteinTrajectory::discretizeTrajectory)
                .def("discretizeTrajectoryH5", &patchyProteinTrajectory::discretizeTrajectoryH5)
                .def("write2H5file", &patchyProteinTrajectory::write2H5file<double, 8>)
                .def("writeChunk2H5file", &patchyProteinTrajectory::writeChunk2H5file<double, 8>);

        // Alternative version of patchyProteinTrajectory
        py::class_<patchyProteinTrajectory2, patchyProteinTrajectory>(m, "patchyProtein2", "alternative discrete "
                                                                "trajectory of patchy protein example(#particles "
                                                                "or #pairs of particles, approx size)")
                .def(py::init<int &, int &>())
                .def(py::init<int &, int &, double &, double &>())
                .def("sampleDiscreteState", &patchyProteinTrajectory2::sampleDiscreteState)
                .def("getState", &patchyProteinTrajectory2::sampleDiscreteState);



        /* Not defined as child class since parent class is a template with virtual fucntions, so need to
         * add all functions from parent classes manually (all the way to the original trajectory.hpp parent).
         * Also note sampleDiscreteState and getState are the same function. */
        py::class_<MAPKtrajectory>(m, "MAPK", "discrete trajectory for MAPK application"
                                              "example(#particles or #pairs of particles, approx size)")
                .def(py::init<int &, int &, double &>())
                .def(py::init<int &, int &, double &, double &, double &>())
                .def(py::init<int &, int &, double &, int &, int &, double &, double &>())
                .def("setAnglePatches", &MAPKtrajectory::setBoundary)
                .def("setBoundary", &MAPKtrajectory::setBoundary)
                .def("setTolerances", &MAPKtrajectory::setTolerances)
                .def("sample", &MAPKtrajectory::sample)
                .def("sampleRelative", &MAPKtrajectory::sampleRelative)
                .def("write2file", &MAPKtrajectory::write2file<double>)
                .def("emptyBuffer", &MAPKtrajectory::emptyBuffer)
                .def("sampleDiscreteTrajectory", &MAPKtrajectory::sampleDiscreteTrajectory)
                .def("sampleDiscreteState", &MAPKtrajectory::sampleDiscreteState)
                .def("getState", &MAPKtrajectory::sampleDiscreteState)
                .def("discretizeTrajectory", &MAPKtrajectory::discretizeTrajectory)
                .def("discretizeTrajectoryH5", &MAPKtrajectory::discretizeTrajectoryH5)
                .def("write2H5file", &MAPKtrajectory::write2H5file<double, 10>)
                .def("writeChunk2H5file", &MAPKtrajectory::writeChunk2H5file<double, 10>);

    }
}
