# CMakeLists.txt per esempio standalone integrazione AstDyn
# Questo file mostra come integrare i moduli FASE 1 in un progetto

cmake_minimum_required(VERSION 3.15)
project(AstDynIntegrationExample VERSION 1.0 LANGUAGES CXX)

# Richiedi C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find dependencies
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(AstDyn REQUIRED)

# Include directories
set(TEMPLATES_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../templates_ioccultcalc")

include_directories(
    ${TEMPLATES_DIR}/include
)

# Source files
set(INTEGRATION_SOURCES
    ${TEMPLATES_DIR}/src/eq1_parser.cpp
    ${TEMPLATES_DIR}/src/orbital_conversions.cpp
    ${TEMPLATES_DIR}/src/astdyn_wrapper.cpp
)

# Create library from integration modules
add_library(astdyn_integration STATIC ${INTEGRATION_SOURCES})

target_include_directories(astdyn_integration
    PUBLIC
        ${TEMPLATES_DIR}/include
)

target_link_libraries(astdyn_integration
    PUBLIC
        AstDyn::AstDyn
        Eigen3::Eigen
)

# Executable: test standalone
add_executable(test_astdyn_integration
    test_astdyn_integration_standalone.cpp
)

target_link_libraries(test_astdyn_integration
    PRIVATE
        astdyn_integration
)

# Compiler flags
target_compile_options(test_astdyn_integration PRIVATE
    -Wall -Wextra -Wpedantic
    $<$<CONFIG:Release>:-O3>
    $<$<CONFIG:Debug>:-g -O0>
)

# Installation
install(TARGETS test_astdyn_integration
    RUNTIME DESTINATION bin
)

# Test
enable_testing()

# Test con asteroid 17030
add_test(
    NAME test_17030_propagation
    COMMAND test_astdyn_integration 
            ${CMAKE_CURRENT_SOURCE_DIR}/../astdyn/data/17030.eq1
            2460643.77083
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

# Summary
message(STATUS "========================================")
message(STATUS "  AstDyn Integration Example")
message(STATUS "========================================")
message(STATUS "Eigen3 version: ${Eigen3_VERSION}")
message(STATUS "AstDyn found: ${AstDyn_FOUND}")
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "========================================")
