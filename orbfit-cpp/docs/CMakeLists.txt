cmake_minimum_required(VERSION 3.15)

# Find Doxygen
find_package(Doxygen)

if(DOXYGEN_FOUND)
    # Set input and output files
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    
    # Copy Doxyfile to build directory
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    
    # Add custom target to build documentation
    add_custom_target(docs
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM
    )
    
    message(STATUS "Doxygen found. Run 'make docs' to generate documentation.")
else()
    message(WARNING "Doxygen not found. Documentation will not be generated.")
endif()
