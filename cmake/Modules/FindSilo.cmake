# - Find Silo, a library for reading and writing unstructured mesh data.
#
# This module will define the following variables:
#  SILO_INCLUDE_DIRS - Location of the SILO includes
#  SILO_LIBRARIES    - Required libraries for the SILO C bindings.
#  SILO_FOUND        - true if SILO was found on the system
#  SILO_LIBRARY_DIRS - the full set of library directories
#------------------------------------------------------------------------
include(SelectLibraryConfigurations)
#include(${CMAKE_CURRENT_LIST_DIR}/FindPackageHandleStandardArgs.cmake)

# Try to find the Silo header file.
find_path( SILO_INCLUDE_DIRS silo.h
           HINTS ${SILO_ROOT} ${HDF5_INCLUDE_DIRS} ${CMAKE_INSTALL_PREFIX} ${possible_inc_dir}
           PATHS ${HDF5_INCLUDE_DIRS}
           PATH_SUFFIXES include)

# Try to find the Silo library file.
find_library( SILO_LIBRARIES
              NAMES siloh5
              HINTS ${SILO_ROOT} ${HDF5_ROOT}/lib ${CMAKE_INSTALL_PREFIX} ${possible_lib_dir}
              PATHS ${HDF5_ROOT}/lib
              PATH_SUFFIXES lib)

if (SILO_INCLUDE_DIRS AND SILO_LIBRARIES)
  set(SILO_FOUND TRUE)
  set(SILO_LIBRARY_DIRS ${SILO_LIBRARIES})
else ()
  set(SILO_FOUND FALSE)
endif ()

find_package_handle_standard_args( Silo
  DEFAULT_MSG
  SILO_LIBRARIES
  SILO_INCLUDE_DIRS)
