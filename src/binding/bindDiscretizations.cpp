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
        /* On binding of spherePartition, change default holder from unique_ptr to shared_ptr to allow msmrdIntegrator to
         * set spherePartition as share pointer. This should also be the case for all of its child classes. */
        py::class_<spherePartition, std::shared_ptr<spherePartition>>(m, "spherePartition", "Equal area spherical "
                                                                                            "partition (numSections)")
                .def(py::init<int &>())
                .def_property_readonly("numSections", &spherePartition::getNumSections)
                .def("getPartition", &spherePartition::getPartition)
                .def("getSectionNumber", &spherePartition::getSectionNumberPyBind)
                .def("getAngles", &spherePartition::getAngles)
                .def("setThetasOffset", &spherePartition::setThetasOffset);

        py::class_<halfSpherePartition, spherePartition, std::shared_ptr<halfSpherePartition>>(m, "halfSpherePartition",
                                                                                   "Equal area partition of the"
                                                                                   "half sphere (numSections)")
                .def(py::init<int &>());


        py::class_<quaternionPartition>(m, "quaternionPartition", "Volumetric partition of unit sphere to discretize "
                                                                  "quaternions (numRadialSecs, numSphericalSecs)")
            .def(py::init<int &, int &>())
            .def_property_readonly("numSections", &quaternionPartition::getNumSections)
            .def("getPartition", &quaternionPartition::getPartition)
            .def("getSectionNumber", &quaternionPartition::getSectionNumberPyBind)
            .def("getSectionIntervals", &quaternionPartition::getSectionIntervals)
            .def("setThetasOffset", &quaternionPartition::setThetasOffset);

        /* On binding of positionOrientationPartition, change default holder from unique_ptr to shared_ptr
         * to allow msmrdIntegrator to set positionOrientationPartition as shared pointer. This should also
         * be done for all of its child classes. */
        py::class_<positionOrientationPartition, std::shared_ptr<positionOrientationPartition>>(m,
                                                    "positionOrientationPartition", "Combines quaternion partion"
                                                    "with spherical partition to discretize relative position and "
                                                    "orientation quaternions (elativeDistanceCutOff, "
                                                    "numSphericalSectionsPos,numRadialSectionsQuat, "
                                                    "numSphericalSectionsQuat)")
                .def(py::init<double &, int &, int &, int &>())
                .def_property_readonly("numSections", &positionOrientationPartition::getNumSections)
                .def("getSectionNumber", &positionOrientationPartition::getSectionNumberPyBind)
                .def("getSectionNumbers", &positionOrientationPartition::getSectionNumbers)
                .def("setThetasOffset", &positionOrientationPartition::setThetasOffset);

    };
}
