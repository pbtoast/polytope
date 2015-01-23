### - Find Silo, a library for reading and writing unstructured mesh data.
###
### This module will define the following variables:
###  SILO_INCLUDE_DIRS - Location of the hdf5 includes
###  SILO_LIBRARIES - Required libraries for the HDF5 C bindings.
###  SILO_FOUND - true if HDF5 was found on the system
###  SILO_LIBRARY_DIRS - the full set of library directories
##
##include(SelectLibraryConfigurations)
##include(${CMAKE_CURRENT_LIST_DIR}/FindPackageHandleStandardArgs.cmake)
##
### Try to find the silodiff executable.
##find_program( silodiff_exe 
##    NAMES silodiff
##    PATH_SUFFIXES bin Bin 
##    DOC "Silodiff program." )
##mark_as_advanced( silodiff_exe )
##
### Intuit directories for silo from silodiff.
##if (silodiff_exe)
##  string(REPLACE "/bin/silodiff" "/include" possible_inc_dir ${silodiff_exe})
##  string(REPLACE "/bin/silodiff" "/lib" possible_lib_dir ${silodiff_exe})
##  message("-- found silodiff at ${silodiff_exe}")
##endif ()
##
### Try to find the Silo header file.
##find_path( SILO_INCLUDE_DIRS silo.h
##           HINTS ${HDF5_INCLUDE_DIRS} ${possible_inc_dir}
##           PATHS ${HDF5_INCLUDE_DIRS}
##           PATH_SUFFIXES include)
##
### Try to find the Silo library file.
##find_library( SILO_LIBRARIES siloh5
##              HINTS ${HDF5_LIBRARY_DIRS} ${possible_lib_dir}
##              PATHS ${HDF5_LIBRARY_DIRS}
##              PATH_SUFFIXES lib)
##if (SILO_INCLUDE_DIRS AND SILO_LIBRARIES)
##  set(SILO_FOUND TRUE)
##else ()
##  set(SILO_FOUND FALSE)
##endif ()
##
##find_package_handle_standard_args( Silo DEFAULT_MSG 
##    SILO_LIBRARIES 
##    SILO_INCLUDE_DIRS
##)
IF(SILO_INCLUDE_DIR)
  IF(SILO_LIBRARY)
    SET(SILO_LIBRARIES ${SILO_LIBRARY})
    SET(SILO_FOUND true)

    #The following deprecated settings are for backwards compatibility with CMake1.4
    SET (SILO_INCLUDE_PATH ${SILO_INCLUDE_DIR})

  ELSE()
    IF(SILO_FIND_REQURIED)
      message(SEND_ERROR "Unable to find the requested SILO libraries.")
    ENDIF()
  ENDIF()
ENDIF()

# handle the QUIETLY and REQUIRED arguments and set SILO_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(SILO 
                                  DEFAULT_MSG 
                                  SILO_LIBRARY 
                                  SILO_INCLUDE_DIR)

MARK_AS_ADVANCED(SILO_INCLUDE_DIR
                 SILO_LIBRARY)
