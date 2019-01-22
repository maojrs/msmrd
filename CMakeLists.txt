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

# Find and add headers location of HDF5 library
# Uncomment to manually specify HDF5 installation directory (useful when multiple installations available)
#set(ENV{HDF5_DIR} "/srv/public/maojrs/miniconda3/envs/msmrd-dev/include/") 
find_package(HDF5 COMPONENTS CXX)
message(STATUS "HDF5_FOUND: ${HDF5_FOUND}, version ${HDF5_VERSION}, directory ${HDF5_INCLUDE_DIRS}")
include_directories(${HDF5_INCLUDE_DIR})

# Set include and source
include_directories(include)
add_subdirectory(libraries/pybind11)
set(bindings_python_version 3.6)
set(SOURCES
        src/msm.cpp
        src/randomgen.cpp
        src/simulation.cpp
        src/boundaries/boundary.cpp
        src/boundaries/box.cpp
        src/boundaries/sphere.cpp
        src/integrators/integrator.cpp
        src/integrators/overdampedLangevin.cpp
        src/integrators/overdampedLangevinMarkovSwitch.cpp
        src/potentials/dipole.cpp
        src/potentials/gaussians3D.cpp
        src/potentials/gayBerne.cpp
        src/potentials/harmonicRepulsion.cpp
        src/potentials/patchyParticle.cpp
        src/potentials/patchyParticleAngular.cpp
        src/potentials/patchyProtein.cpp
        src/potentials/patchyProteinMarkovSwitch.cpp
        src/potentials/potentials.cpp
        src/trajectories/trajectory.cpp
        )

set(PY_SOURCES
        src/binding/binding.hpp
        src/binding/binding.cpp
        src/binding/bindBoundaries.cpp
        src/binding/bindIntegrators.cpp
        src/binding/bindInternal.cpp
        src/binding/bindMarkovModels.cpp
        src/binding/bindParticles.cpp
        src/binding/bindPotentials.cpp
        src/binding/bindSimulation.cpp
        src/binding/bindTrajectory.cpp
        include/msm.hpp
        include/particle.hpp
        include/quaternion.hpp
        include/randomgen.hpp
        include/simulation.hpp
        include/tools.hpp
        include/vec3.hpp
        include/boundaries/boundary.hpp
        include/boundaries/box.hpp
        include/boundaries/sphere.hpp
        include/integrators/integrator.hpp
        include/integrators/overdampedLangevin.hpp
        include/integrators/overdampedLangevinMarkovSwitch.hpp
        include/potentials/potentials.hpp
        include/potentials/dipole.hpp
        include/potentials/gaussians3D.hpp
        include/potentials/gayBerne.hpp
        include/potentials/harmonicRepulsion.hpp
        include/potentials/patchyParticle.hpp
        include/potentials/patchyParticleAngular.hpp
        include/potentials/patchyProtein.hpp
        include/potentials/patchyProteinMarkovSwitch.hpp
        include/trajectories/trajectory.hpp
        )

add_library(msmrd2core SHARED ${SOURCES})

target_link_libraries(msmrd2core ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})

#target_include_directories(msmrd2core PUBLIC include libraries/pybind11/include)
pybind11_add_module(msmrd2binding MODULE ${PY_SOURCES})

#target_include_directories(msmrd2 PUBLIC include libraries/pybind11/include)
target_link_libraries(msmrd2binding PUBLIC msmrd2core)


#EXECUTE_PROCESS(COMMAND git submodule update --init WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# Add catch2 library
add_library(catch2 INTERFACE)
target_include_directories(catch2 INTERFACE "libraries/Catch2/single_include")
add_subdirectory(tests/cpp)
