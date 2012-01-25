# Minimum CMake version.
cmake_minimum_required (VERSION 2.8.5)

# Adjust Polytope's module path.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Set compilers. This must be done before enabling languages.
set(CMAKE_C_COMPILER "${CC}")
message("-- C compiler is ${CMAKE_C_COMPILER}")
set(CMAKE_CXX_COMPILER "${CXX}")
message("-- C++ compiler is ${CMAKE_CXX_COMPILER}")

# Build everything as static libs.
set (BUILD_SHARED_LIBS OFF)

# Project and version numbers.
project (polytope)
set (POLYTOPE_VERSION_MAJOR 1)
set (POLYTOPE_VERSION_MINOR 0)

# Figure out the system type.
if (${APPLE})
  set(SYS_FLAGS "-DAPPLE=1")
else ()
  if (${LINUX})
    set(SYS_FLAGS "-DLINUX=1")
  endif ()
endif ()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SYS_FLAGS}")

# General compiler flags.
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-sign-compare -ansi")
endif ()

# Figure out MPI.
if (USE_MPI EQUAL 1)
  # CC and CXX should already have been set in Makefile or wherever.
  set(MPIEXEC mpirun)
  set(MPIEXEC_NUMPROC_FLAG -np)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_MPI")
endif ()

# Find HDF5.
set(HDF5_USE_STATIC_LIBRARIES TRUE)
find_package(HDF5)
if (${HDF5_FOUND})
  # Do we need parallel HDF5?
  if (USE_MPI EQUAL 1)
    if (NOT ${HDF5_IS_PARALLEL})
      message("-- Parallel HDF5 not found. Skipping.")
      set(HDF5_FOUND FALSE)
    endif ()
  endif ()
endif ()
if (${HDF5_FOUND})
  include_directories("${HDF5_INCLUDE_DIRS}")
  link_directories("${HDF5_LIBRARY_DIRS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_HDF5")

  # Get rid of 'debug' and 'optimized' detritis.
  list(REMOVE_ITEM HDF5_LIBRARIES debug optimized)
endif ()

# If HDF5 is here, look for Silo.
if (${HDF5_FOUND})
  find_package(Silo)
  if (${SILO_FOUND})
    include_directories("${SILO_INCLUDE_DIRS}")
#   FIXME: For now, we assume libsiloh5 lives with libhdf5.
#    link_directories("${SILO_LIBRARY_DIRS}")
    set(CHARYBDIS_TEST_LIBS "${CHARYBDIS_TEST_LIBS} ${SILO_LIBRARIES} ${HDF5_LIBRARIES}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_SILO")
  else ()
    set(CHARYBDIS_TEST_LIBS "${CHARYBDIS_TEST_LIBS} ${HDF5_LIBRARIES}")
  endif ()
else ()
  set(SILO_FOUND FALSE)
endif ()

# Build a configuration header file from our options.
configure_file(
  "${PROJECT_SOURCE_DIR}/src/polytope.hh.in"
  "${PROJECT_BINARY_DIR}/polytope.hh"
)

# Turn on testing support.
enable_testing()

# Include the binary directory in the header file search path.
include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/src")

# Libraries

# Voro 2d/3d libraries.
include_directories(voro_2d)
add_subdirectory(voro_2d)
include_directories(voro_3d)
add_subdirectory(voro_3d)

# Polytope proper.
include_directories("${PROJECT_SOURCE_DIR}/src/")
add_subdirectory(src)

# Tests directory
add_subdirectory(tests)

# Library install targets
install (FILES "${PROJECT_BINARY_DIR}/polytope.hh"        
         DESTINATION include)

# build a CPack driven installer package
include(InstallRequiredSystemLibraries)
set (CPACK_RESOURCE_FILE_LICENSE  
     "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set (CPACK_PACKAGE_VERSION_MAJOR "${POLYTOPE_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${POLYTOPE_VERSION_MINOR}")
include(CPack)

# enable dashboard scripting
include (CTest)
  