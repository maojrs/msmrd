cmake_minimum_required(VERSION 2.8.12)
project(msmrd2)

#set c++14 std for old compilers
macro(use_cxx14)
    if (CMAKE_VERSION VERSION_LESS "6.3")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
    endif ()
endmacro(use_cxx14)
use_cxx14()

FIND_PROGRAM(PYTHON "python")

include_directories(include)
add_subdirectory(libraries/pybind11)
set(bindings_python_version 3.6)
set(SOURCES
        src/boundary.cpp
        src/msm.cpp
        src/randomgen.cpp
        src/simulation.cpp
        src/integrators/integrator.cpp
        src/integrators/odLangevin.cpp
        src/integrators/odLangevinMarkovSwitch.cpp
        src/potentials/potentials.cpp
        src/potentials/dipole.cpp
        src/potentials/gaussians3D.cpp
        src/potentials/harmonicRepulsion.cpp
        src/potentials/gayBerne.cpp
)

set(PY_SOURCES
        src/binding/binding.hpp
        src/binding/binding.cpp
        src/binding/bindBoundaries.cpp
        src/binding/bindIntegrators.cpp
        src/binding/bindMarkovModels.cpp
        src/binding/bindParticles.cpp
        src/binding/bindPotentials.cpp
        include/boundary.hpp
        include/msm.hpp
        include/particle.hpp
        include/quaternion.hpp
        include/randomgen.hpp
        include/simulation.hpp
        include/vec3.hpp
        include/integrators/integrator.hpp
        include/integrators/odLangevin.hpp
        include/integrators/odLangevinMarkovSwitch.hpp
        include/potentials/potentials.hpp
        include/potentials/dipole.hpp
        include/potentials/gaussians3D.hpp
        include/potentials/harmonicRepulsion.hpp
        include/potentials/gayBerne.hpp
        )

add_library(msmrd2core SHARED ${SOURCES})
#target_include_directories(msmrd2core PUBLIC include libraries/pybind11/include)

pybind11_add_module(msmrd2binding MODULE ${PY_SOURCES})

#target_include_directories(msmrd2 PUBLIC include libraries/pybind11/include)
target_link_libraries(msmrd2binding PUBLIC msmrd2core)



#EXECUTE_PROCESS(COMMAND git submodule update --init WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_library(catch2 INTERFACE)
target_include_directories(catch2 INTERFACE "libraries/Catch2/single_include")
add_subdirectory(tests/cpp)
