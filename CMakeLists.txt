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
        src/simulation.cpp
        src/integrators/integrator.cpp
        src/msm.cpp
        src/randomgen.cpp
        src/integrators/integrator.cpp
        src/integrators/odLangevin.cpp
        src/integrators/odLangevinMarkovSwitch.cpp
        src/potentials/gaussians3D.cpp
        src/potentials/harmonicRepulsion.cpp
        src/potentials/gayBerne.cpp
        )

set(PY_SOURCES
        src/binding/binding.cpp
        include/quaternion.hpp
        include/vec3.hpp
        include/particle.hpp
        include/simulation.hpp
        include/msm.hpp
        include/randomgen.hpp
        )

set(PY_SOURCES_INTEGRATORS
        src/binding/bindIntegrators.cpp
        include/integrators/integrator.hpp
        include/integrators/odLangevin.hpp
        include/integrators/odLangevinMarkovSwitch.hpp
        )

set(PY_SOURCES_POTENTIALS
        src/binding/bindPotentials.cpp
        include/potentials/potential.hpp
        include/potentials/gaussians3D.hpp
        include/potentials/harmonicRepulsion.hpp
        include/potentials/gayBerne.hpp
        )

add_library(msmrd2core SHARED ${SOURCES})
#target_include_directories(msmrd2core PUBLIC include libraries/pybind11/include)

pybind11_add_module(main MODULE ${PY_SOURCES})
pybind11_add_module(integrators MODULE ${PY_SOURCES_INTEGRATORS})
pybind11_add_module(potentials MODULE ${PY_SOURCES_POTENTIALS})


#target_include_directories(msmrd2 PUBLIC include libraries/pybind11/include)
target_link_libraries(main PUBLIC msmrd2core)
target_link_libraries(integrators PUBLIC msmrd2core)
target_link_libraries(potentials PUBLIC msmrd2core)



#EXECUTE_PROCESS(COMMAND git submodule update --init WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_library(catch2 INTERFACE)
target_include_directories(catch2 INTERFACE "libraries/Catch2/single_include")
add_subdirectory(tests/cpp)
