cmake_minimum_required(VERSION 3.15)

# Test executable for Phase 1
add_executable(orbfit_tests
    test_constants.cpp
    test_types.cpp
)

target_link_libraries(orbfit_tests
    PRIVATE
        orbfit
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(orbfit_tests)

# Test executable for Phase 2 - Math
add_executable(orbfit_math_tests
    test_math.cpp
)

target_link_libraries(orbfit_math_tests
    PRIVATE
        orbfit
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(orbfit_math_tests)

# Test executable for Phase 2 - Time
add_executable(orbfit_time_tests
    test_time.cpp
)

target_link_libraries(orbfit_time_tests
    PRIVATE
        orbfit
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(orbfit_time_tests)

# Test executable for Phase 2 - Utils
add_executable(orbfit_utils_tests
    test_utils.cpp
)

target_link_libraries(orbfit_utils_tests
    PRIVATE
        orbfit
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(orbfit_utils_tests)

# Test executable for Phase 3 - Coordinates
add_executable(orbfit_coordinates_tests
    test_coordinates.cpp
)

target_link_libraries(orbfit_coordinates_tests
    PRIVATE
        orbfit
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(orbfit_coordinates_tests)

# Test executable for Phase 3 - Reference Frames
add_executable(orbfit_reference_frames_tests
    test_reference_frames.cpp
)

target_link_libraries(orbfit_reference_frames_tests
    PRIVATE
        orbfit
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(orbfit_reference_frames_tests)

# Test executable for Phase 3 - Jacobian Matrices
add_executable(orbfit_jacobian_tests
    test_jacobians.cpp
)

target_link_libraries(orbfit_jacobian_tests
    PRIVATE
        orbfit
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(orbfit_jacobian_tests)

# Test executable for Phase 6 - Propagation
add_executable(orbfit_propagation_tests
    test_propagation.cpp
)

target_link_libraries(orbfit_propagation_tests
    PRIVATE
        orbfit
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(orbfit_propagation_tests)

# Test executable for Phase 3 - Validation Against Known Data
add_executable(orbfit_validation_tests
    test_validation.cpp
)

target_link_libraries(orbfit_validation_tests
    PRIVATE
        orbfit
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(orbfit_validation_tests)

# ============================================================================
# Ephemeris Tests
# ============================================================================

add_executable(orbfit_ephemeris_tests
    test_ephemeris.cpp
)

target_link_libraries(orbfit_ephemeris_tests
    PRIVATE
        orbfit
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(orbfit_ephemeris_tests)

# ============================================================================
# RKF78 Debug Tests
# ============================================================================

add_executable(orbfit_rkf78_debug
    test_rkf78_debug.cpp
)

target_link_libraries(orbfit_rkf78_debug
    PRIVATE
        orbfit
        GTest::gtest
        GTest::gtest_main
)

gtest_discover_tests(orbfit_rkf78_debug)
