# Find qhull header and library.
#

# This module defines the following uncached variables:
#  QHULL_FOUND, if false, do not try to use qhull.
#  QHULL_INCLUDE_DIRS, where to find the Qhull.h C++ interface header.
#  QHULL_LIBRARIES, the libraries to link against to use the qhull library
#  QHULL_LIBRARY_DIRS, the directory where the qhull library is found.

find_path(
  QHULL_INCLUDE_DIR Qhull.h
  /usr/local/include
  /usr/include
)

if( QHULL_INCLUDE_DIR )
  find_library(
    QHULL_LIBRARY
    NAMES qhull
    PATHS /usr/local/lib /usr/lib
  )
  find_library(
    QHULLCPP_LIBRARY
    NAMES qhullcpp
    PATHS /usr/local/lib /usr/lib
  )
  if( QHULL_LIBRARY AND QHULLCPP_LIBRARY)
    set(QHULL_LIBRARY_DIR "")
    get_filename_component(QHULL_LIBRARY_DIRS ${QHULL_LIBRARY} PATH)
    # Set uncached variables as per standard.
    set(QHULL_FOUND ON)
    set(QHULL_INCLUDE_DIRS ${QHULL_INCLUDE_DIR})
    set(QHULL_LIBRARIES ${QHULLCPP_LIBRARY};${QHULL_LIBRARY})
  endif()
endif()
	    
if(QHULL_FOUND)
  if(NOT QHULL_FIND_QUIETLY)
    message(STATUS "FindQHull: Found Qhull.h and libraries")
  endif(NOT QHULL_FIND_QUIETLY)
else(QHULL_FOUND)
  if(QHULL_FIND_REQUIRED)
    message(FATAL_ERROR "FindQHull: Could not find Qhull.h and/or libraries")
  endif(QHULL_FIND_REQUIRED)
endif(QHULL_FOUND)

