//
// Created by maojrs on 4/2/19.
//

#include "binding.hpp"
#include "discretizations/halfSpherePartition.hpp"
#include "discretizations/quaternionPartition.hpp"
#include "discretizations/spherePartition.hpp"
#include "discretizations/positionOrientationPartition.hpp"

namespace msmrd {
    /*
     * pyBinders for the c++ discretizations classes
     */
    void bindDiscretizations(py::module &m) {
        py::class_<spherePartition>(m, "spherePartition", "Equal area spherical partition (numSections)")
                .def(py::init<int &>())
                .def_property_readonly("numSections", &spherePartition::getNumSections)
                .def("getPartition", &spherePartition::getPartition)
                .def("getSectionNumber", &spherePartition::getSectionNumberPyBind)
                .def("getAngles", &spherePartition::getAngles);

        py::class_<halfSpherePartition, spherePartition>(m, "halfSpherePartition", "Equal area partition of the"
                                                                                   "half sphere (numSections)")
                .def(py::init<int &>());


        py::class_<quaternionPartition>(m, "quaternionPartition", "Volumetric partition of unit sphere to discretize "
                                                                  "quaternions (numRadialSecs, numSphericalSecs)")
            .def(py::init<int &, int &>())
            .def_property_readonly("numSections", &quaternionPartition::getNumSections)
            .def("getPartition", &quaternionPartition::getPartition)
            .def("getSectionNumber", &quaternionPartition::getSectionNumberPyBind)
            .def("getSectionIntervals", &quaternionPartition::getSectionIntervals);

        py::class_<positionOrientationPartition>(m, "positionOrientationPartition", "Combines quaternion partion"
                                                    "with spherical partition to discretize relative position and "
                                                    "orientation quaternions (elativeDistanceCutOff, "
                                                    "numSphericalSectionsPos,numRadialSectionsQuat, "
                                                    "numSphericalSectionsQuat)")
                .def(py::init<double &, int &, int &, int &>())
                .def_property_readonly("numSections", &positionOrientationPartition::getNumSections)
                .def("getSectionNumber", &positionOrientationPartition::getSectionNumberPyBind)
                .def("getSection", &positionOrientationPartition::getSection);

    };
}
