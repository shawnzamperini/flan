#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "flan" for configuration ""
set_property(TARGET flan APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(flan PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "/home/zamp/github/flan/lib/libflan.so"
  IMPORTED_SONAME_NOCONFIG "libflan.so"
  )

list(APPEND _cmake_import_check_targets flan )
list(APPEND _cmake_import_check_files_for_flan "/home/zamp/github/flan/lib/libflan.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
