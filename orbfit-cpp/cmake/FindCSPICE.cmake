# Try to find CSPICE library
# Once done this will define:
#  CSPICE_FOUND - System has CSPICE
#  CSPICE_INCLUDE_DIRS - The CSPICE include directories
#  CSPICE_LIBRARIES - The libraries needed to use CSPICE

find_path(CSPICE_INCLUDE_DIR 
    NAMES SpiceUsr.h
    PATHS
        /usr/local/include/cspice
        /usr/include/cspice
        $ENV{CSPICE_ROOT}/include
        ${CSPICE_ROOT}/include
)

find_library(CSPICE_LIBRARY
    NAMES cspice
    PATHS
        /usr/local/lib
        /usr/lib
        $ENV{CSPICE_ROOT}/lib
        ${CSPICE_ROOT}/lib
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CSPICE DEFAULT_MSG
    CSPICE_LIBRARY CSPICE_INCLUDE_DIR
)

if(CSPICE_FOUND)
    set(CSPICE_LIBRARIES ${CSPICE_LIBRARY})
    set(CSPICE_INCLUDE_DIRS ${CSPICE_INCLUDE_DIR})
    
    if(NOT TARGET CSPICE::CSPICE)
        add_library(CSPICE::CSPICE UNKNOWN IMPORTED)
        set_target_properties(CSPICE::CSPICE PROPERTIES
            IMPORTED_LOCATION "${CSPICE_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${CSPICE_INCLUDE_DIR}"
        )
    endif()
endif()

mark_as_advanced(CSPICE_INCLUDE_DIR CSPICE_LIBRARY)
