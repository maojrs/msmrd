//
// Created by maojrs on 9/7/18.
//
#include "binding.hpp"
#include "markovModels/discreteTimeMarkovModel.hpp"
#include "markovModels/continuousTimeMarkovModel.hpp"
#include "markovModels/msmrdMarkovModel.hpp"

namespace msmrd {
    // Aliases for classes with long names.
    using msm = msmrd::discreteTimeMarkovStateModel;
    using ctmsm = msmrd::continuousTimeMarkovStateModel;
    using msmrdMSM = msmrd::msmrdMarkovStateModel;
    /*
     * pyBinders for the c++ Markov state models (MSMs) classes
     */
    void bindMarkovModels(py::module &m) {
        py::class_<msm, markovModel>(m, "discreteTimeMarkovStateModel", "discrete time Markov state model (MSM ID, "
                                                           "transition matrix, lagtime, seed)")
                .def(py::init<int &, std::vector<std::vector<double>> &, double &, long &>())
                .def("propagate", &msm::propagate)
                .def("propagateMSM", &msm::propagateMSM);


        py::class_<ctmsm, markovModel>(m, "continuousTimeMarkovStateModel", "continuous time Markov state model (MSM ID, "
                                                               "transition rate matrix, seed)")
                .def(py::init<int &, std::vector<std::vector<double>> &, long &>())
                .def("propagate", &ctmsm::propagate)
                .def("propagateMSM", &ctmsm::propagateMSM);

        py::class_<msmrdMSM>(m, "msmrdMarkovStateModel", "continuous time Markov state model specialized to use with"
                                                            "MSM/RD integration (num. of bound states, num. transition"
                                                            "states, seed, rate dictionary")
                .def(py::init<unsigned int &, unsigned int &, long &, std::map<std::string, float> &>())
                .def("getRate", &msmrdMSM::getRate)
                .def("computeTransition2BoundState", &msmrdMSM::computeTransition2BoundState)
                .def("computeTransitionFromoundState", &msmrdMSM::computeTransitionFromBoundState)
                .def("setDbound", &msmrdMSM::setDbound)
                .def("setMaxNumberBoundStates", &msmrdMSM::setMaxNumberBoundStates);
    }

}
