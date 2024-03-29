//
// Created by maojrs on 9/7/18.
//
#include <integrators/msmrdIntegrator.hpp>
#include "binding.hpp"
#include "markovModels/discreteTimeMarkovModel.hpp"
#include "markovModels/continuousTimeMarkovModel.hpp"
#include "markovModels/msmrdMarkovModel.hpp"

namespace msmrd {
    // Aliases for classes with long names.
    using msm = msmrd::discreteTimeMarkovStateModel;
    using ctmsm = msmrd::continuousTimeMarkovStateModel;
    using msmrdMSM = msmrd::msmrdMarkovModel;

    /*
     * pyBinders for the c++ Markov state models (MSMs) classes
     */
    void bindMarkovModels(py::module &m) {
        py::class_<msm, markovModel>(m, "discreteTimeMarkovStateModel", "discrete time Markov state model (MSM ID, "
                                                           "transition matrix, lagtime, seed)")
                .def(py::init<int &, double &, long &>())
                .def(py::init<int &, std::vector<std::vector<double>> &, double &, long &>())
                .def("propagate", &msm::propagate)
                .def("propagateMSM", &msm::propagateMSM);


        py::class_<ctmsm, markovModel>(m, "continuousTimeMarkovStateModel", "continuous time Markov state model (MSM ID, "
                                                               "transition rate matrix, seed)")
                .def(py::init<int &, long &>())
                .def(py::init<int &, std::vector<std::vector<double>> &, long &>())
                .def("propagate", &ctmsm::propagate)
                .def("propagateMSM", &ctmsm::propagateMSM);


        py::class_<msmrdMSM, msm>(m, "msmrdMarkovModel", "discrete time Markov state model "
                                        "specialized to  use with MSM/RD integration (num. of bound states, "
                                        "max number of bound states, transition probability matrix, MSM active set,"
                                        "seed")
                .def(py::init<unsigned int &, unsigned int &, std::vector<std::vector<double>> &,
                        std::vector<int> &, float &, long &>())
                .def("calculateTransition", &msmrdMSM::calculateTransition)
                .def("setDbound", &msmrdMSM::setDbound)
                .def("setMaxNumberBoundStates", &msmrdMSM::setMaxNumberBoundStates);

    }

}
