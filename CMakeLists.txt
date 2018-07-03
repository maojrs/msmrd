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

#include_directories(include)
include_directories(include)
add_subdirectory(libraries/pybind11)
set(bindings_python_version 3.6)
set(SOURCES
        src/simulation.cpp
        src/integrator.cpp
#        include/quaternion.hpp
#        include/vec3.hpp
#        include/particle.hpp
#        include/simulation.hpp
#        include/integrator.hpp
        )

set(PY_SOURCES
        src/binding.cpp
        include/quaternion.hpp
        include/vec3.hpp
        include/particle.hpp
        include/simulation.hpp
        include/integrator.hpp
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
