#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "ITALOccultLibrary::italoccultlib" for configuration "Release"
set_property(TARGET ITALOccultLibrary::italoccultlib APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(ITALOccultLibrary::italoccultlib PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libitaloccultlib.a"
  )

list(APPEND _cmake_import_check_targets ITALOccultLibrary::italoccultlib )
list(APPEND _cmake_import_check_files_for_ITALOccultLibrary::italoccultlib "${_IMPORT_PREFIX}/lib/libitaloccultlib.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
