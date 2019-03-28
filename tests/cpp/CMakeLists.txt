set(TEST_SOURCES
        testMain.cpp
        testDiscretizations.cpp
        testMarkovModels.cpp
        testPotentials.cpp
        testTrajectories.cpp
        )

add_executable(test_executable ${TEST_SOURCES})

target_link_libraries(test_executable catch2 msmrd2core)

target_include_directories(test_executable PUBLIC
        include/quaternion.h
        include/vec3.h
        include/particle.h
        include/simulation.h
        include/integrator.h)
