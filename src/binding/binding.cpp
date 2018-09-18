
#include "binding.hpp"

/**
 * Bind msmrd2 modules and submodules (functions defined in binding.hpp
 * and implemented in src/bind****.py)
 */
PYBIND11_MODULE(msmrd2binding, module) {

    module.attr("__name__") = "msmrd2";
    module.doc() =  "msmrd with python bindings";

    // Load classes in main module
    msmrd::bindBoundaries(module);
    msmrd::bindMarkovModels(module);
    msmrd::bindParticles(module);

    // Load submodules (integrators and potentials)
    auto integratorsSubmodule = module.def_submodule("integrators");
    msmrd::bindIntegrators(integratorsSubmodule);
    auto potentialsSubmodule = module.def_submodule("potentials");
    msmrd::bindPotentials(potentialsSubmodule);

}
