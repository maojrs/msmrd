set(TEST_SOURCES
        testMain.cpp
        testBoundaries.cpp
        testDiscretizations.cpp
        testIntegrators.cpp
        testMarkovModels.cpp
        testPotentials.cpp
        testTools.cpp
        testTrajectories.cpp)

add_executable(test_executable ${TEST_SOURCES})

target_link_libraries(test_executable catch2 msmrd2core)

target_include_directories(test_executable PUBLIC
        include/quaternion.hpp
        include/vec3.hpp
        include/particle.hpp
        include/simulation.hpp
        include/integrator.hpp)
