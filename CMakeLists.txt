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
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-sign-compare -ansi -fPIC")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-sign-compare")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "XL")
  set(USING_XLC true)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qPIC")
else ()
  set(USING_XLC false)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
endif ()

# Figure out MPI.
if (USE_MPI EQUAL 1)
  # CC and CXX should already have been set in Makefile or wherever.
  set(HAVE_MPI true)
  find_program(SLURM_SRUN_COMMAND srun 
    DOC "Path to Slurm's srun executable")
  if(SLURM_SRUN_COMMAND)
    set(MPIEXEC ${SLURM_SRUN_COMMAND})
    set(MPIEXEC_NUMPROC_FLAG -n)
    set(HAVE_MPIEXEC true)
    message("-- Slurm srun detected. Using ${MPIEXEC}")
  else()
    find_program(MPIRUN_COMMAND mpirun)
    if(MPIRUN_COMMAND)
      set(MPIEXEC mpirun)
      set(MPIEXEC_NUMPROC_FLAG -np)
      set(HAVE_MPIEXEC true)
    else()
      set(MPIEXEC "")
      set(MPIEXEC_NUMPROC_FLAG "")
      set(HAVE_MPIEXEC false)
    endif()
  endif()
else ()
  set(HAVE_MPI false)
  set(HAVE_MPIEXEC false)
endif ()


# Find Boost.
# set(Boost_DEBUG 1)
set(Boost_ADDITIONAL_VERSIONS "1.51" "1.51.0" "1.52" "1.52.0" "1.53" "1.53.0")
find_package(Boost "1.50")
if (USING_XLC) 
  set(Boost_FOUND false)
  set(HAVE_BOOST false)
  message("-- Boost being disabled on IBM for now due to incompatibility with xlC")
endif()
if (Boost_FOUND)
  include_directories("${Boost_INCLUDE_DIR}")
  set(HAVE_BOOST true)
  message("-- Found BOOST installation at ${Boost_INCLUDE_DIR}")
  if (Boost_MINOR_VERSION GREATER 51 AND
      EXISTS ${Boost_INCLUDE_DIR}/boost/polygon/voronoi.hpp)
    set(HAVE_BOOST_VORONOI true)
    message("-- Found BOOST VORONOI in ${Boost_INCLUDE_DIR}/boost/polygon")
  else ()
    set(HAVE_BOOST_VORONOI false)
    message("-- Boost Voronoi not found. BoostTessellator is disabled.")
  endif ()
else ()
  message("-- Unable to find required package BOOST.")
  set(HAVE_BOOST false)
  set(HAVE_BOOST_VORONOI false)
endif ()

if (Boost_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBoost_FOUND=1")
endif()

# Find Python
if (USE_PYTHON EQUAL 1)
  # Do you have the python executable and version number?
  if (PYTHON_EXE AND PYTHON_VERSION)
    string(REPLACE "/bin/python" "" PYTHON_ROOT ${PYTHON_EXE})
    set(PYTHON_INCLUDE_DIR "${PYTHON_ROOT}/include/python${PYTHON_VERSION}")
    set(PYTHON_LIB_DIR "${PYTHON_ROOT}/lib")
    # Check that Pybindgen is installed there
    if (EXISTS ${PYTHON_LIB_DIR}/python${PYTHON_VERSION}/site-packages/pybindgen
	OR     ${PYTHON_LIB_DIR}/python/site-packages/pybindgen)
      # Check that you have Python.h and libpython<version>.a
      if (EXISTS "${PYTHON_INCLUDE_DIR}/Python.h" AND
	  EXISTS "${PYTHON_LIB_DIR}/libpython${PYTHON_VERSION}.a")
	set(PYBINDGEN_DIR "${PROJECT_SOURCE_DIR}/src/PBGWraps")
	set(PYBINDGEN_MODULE_NAME "PolytopeModules")
	set(HAVE_PYTHON true)
	message("-- Found Python version ${PYTHON_VERSION}")
	message("-- Using Python executable ${PYTHON_EXE}")
      else()
	set(HAVE_PYTHON false)
	message("-- Failed to find optional package Python. Bindings will not be generated.")
      endif()
    else()
      set(HAVE_PYTHON false)
      message("-- Failed to find PyBindGen module in ${PYTHON_LIB_DIR}/python${PYTHON_VERSION}/site-packages/")
      message("-- Failed to find PyBindGen module in ${PYTHON_LIB_DIR}/python/site-packages/")
    endif()
  else()
    set(HAVE_PYTHON false)
    message("-- Python executable and vesion number not specified")
  endif()
endif()

#Check for HDF5 and Silo
if (${USE_SILO} EQUAL 1)
  # Find HDF5.
  set(HDF5_USE_STATIC_LIBRARIES true)
  find_package(HDF5)
  if (${HDF5_FOUND})
    # Get rid of'debug' and 'optimized' detritis.
    list(REMOVE_ITEM HDF5_LIBRARIES debug optimized)
    list(REMOVE_ITEM HDF5_INCLUDE_DIRS debug optimized)
    message("-- Found hdf5 files in ${HDF5_INCLUDE_DIRS}")
    
    # Do we have parallel HDF5?
    if (USE_MPI EQUAL 1)
      if (${HDF5_IS_PARALLEL})
	set(HAVE_PARALLEL_HDF5 true)
      else ()
	set(HAVE_PARALLEL_HDF5 false)
      endif ()
    endif ()
  else()
    message("-- hdf5 not found. Skipping.")
  endif()
  
  # If HDF5 is here, look for SILO
  if (${HDF5_FOUND})
    find_package(Silo)
    if (${SILO_FOUND})
      message("-- Found Silo files ${SILO_INCLUDE_DIRS}")
    else()
      message("-- Silo not found. Skipping.")
    endif()
  endif()

  # Add paths only if both SILO and HDF5 found.
  if (HDF5_FOUND AND SILO_FOUND)
    set(HAVE_SILO true)
    set(HAVE_HDF5 true)
    include_directories("${SILO_INCLUDE_DIRS}")
    include_directories("${HDF5_INCLUDE_DIRS}")
#         FIXME: For now, we assume libsiloh5 lives with libhdf5.
#         link_directories("${SILO_LIBRARY_DIRS}")
    link_directories("${HDF5_LIBRARY_DIRS}")
  else()
    message("-- Did not find both Silo and hdf5 files. Skipping.")
    set(HAVE_SILO false)
    set(HAVE_HDF5 false)
    set(SILO_LIBRARIES "")
    set(HDF5_LIBRARIES "")
  endif()
else()
  set(HAVE_SILO false)
  set(HAVE_HDF5 false)
  set(SILO_LIBRARIES "")
  set(HDF5_LIBRARIES "")
  message("-- Not using SILO and HDF5.")
endif()


# Enable Triangle if triangle.h and triangle.c are found.
# We disable this on IBM xlC since Boost.Geometry does not compile correctly there for now.
set(CWD ${CMAKE_CURRENT_SOURCE_DIR})
if (EXISTS ${CWD}/src/triangle.h AND 
    EXISTS ${CWD}/src/triangle.c AND
    Boost_FOUND)
  set(TRIANGLE_FOUND true)
  set(HAVE_TRIANGLE true)
else()
  set(TRIANGLE_FOUND false)
  set(HAVE_TRIANGLE false)
endif()

# Enable Tetgen if the proper files are found.
if (EXISTS ${CWD}/src/tetgen.h AND 
    EXISTS ${CWD}/src/tetgen.cxx)
  set(TETGEN_FOUND true)
  set(HAVE_TETGEN true)
else()
  set(TETGEN_FOUND false)
  set(HAVE_TETGEN false)
endif()

# Build a configuration header file from our options.
configure_file(
  "${PROJECT_SOURCE_DIR}/src/polytope.hh.in"
  "${PROJECT_BINARY_DIR}/polytope.hh"
)


if (HEADER_ONLY EQUAL 0)
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
  if (${TESTING} EQUAL 1)
    enable_testing()
    add_subdirectory(tests)
    include(CTest)
  endif()
  
  # Polytope C library.
  if (BUILD_C_INTERFACE EQUAL 1)
    include_directories("${PROJECT_SOURCE_DIR}/polytope_c/")
    add_subdirectory(polytope_c)
  endif()

  # Python bindings.
  if(HAVE_PYTHON)
    include_directories("${PROJECT_SOURCE_DIR}/src/PBGWraps")
    include_directories("${PROJECT_SOURCE_DIR}/src/PBGWraps/")
    add_subdirectory(src/PBGWraps)
  endif()
endif()

# Library install targets
install (FILES "${PROJECT_BINARY_DIR}/polytope.hh" 
               "${PROJECT_SOURCE_DIR}/src/Tessellator.hh" 
               "${PROJECT_SOURCE_DIR}/src/TessellatorInline.hh" 
               "${PROJECT_SOURCE_DIR}/src/Tessellation.hh" 
               "${PROJECT_SOURCE_DIR}/src/PLC.hh" 
	       "${PROJECT_SOURCE_DIR}/src/PLC_CSG_2d.hh"
	       "${PROJECT_SOURCE_DIR}/src/PLC_CSG_3d.hh"
               "${PROJECT_SOURCE_DIR}/src/ReducedPLC.hh"
               "${PROJECT_SOURCE_DIR}/src/simplifyPLCfacets.hh"
               "${PROJECT_SOURCE_DIR}/src/QuantTessellation.hh" 
	       "${PROJECT_SOURCE_DIR}/src/QuantizedCoordinates.hh"
               "${PROJECT_SOURCE_DIR}/src/VoroPP_2d.hh" 
               "${PROJECT_SOURCE_DIR}/src/VoroPP_3d.hh" 
               "${PROJECT_SOURCE_DIR}/src/TriangleTessellator.hh"
	       "${PROJECT_SOURCE_DIR}/src/BoostTessellator.hh"
               "${PROJECT_SOURCE_DIR}/src/TetgenTessellator.hh"
	       "${PROJECT_SOURCE_DIR}/src/MeshEditor.hh"
               "${PROJECT_SOURCE_DIR}/src/polytope_internal.hh" 
               "${PROJECT_SOURCE_DIR}/src/polytope_geometric_utilities.hh" 
               "${PROJECT_SOURCE_DIR}/src/polytope_parallel_utilities.hh"
               "${PROJECT_SOURCE_DIR}/src/polytope_serialize.hh" 
	       "${PROJECT_SOURCE_DIR}/src/polytope_tessellator_utilities.hh"
               "${PROJECT_SOURCE_DIR}/src/polytope_write_OOGL.hh" 
	       "${PROJECT_SOURCE_DIR}/src/OrphanageBase.hh"
	       "${PROJECT_SOURCE_DIR}/src/BoostOrphanage.hh"
	       "${PROJECT_SOURCE_DIR}/src/Clipper2d.hh"
               "${PROJECT_SOURCE_DIR}/src/ErrorHandler.hh" 
               "${PROJECT_SOURCE_DIR}/src/convexHull_2d.hh" 
               "${PROJECT_SOURCE_DIR}/src/convexHull_3d.hh" 
               "${PROJECT_SOURCE_DIR}/src/Point.hh" 
               "${PROJECT_SOURCE_DIR}/src/KeyTraits.hh" 
	       "${PROJECT_SOURCE_DIR}/src/DimensionTraits.hh"
               "${PROJECT_SOURCE_DIR}/src/SiloWriter.hh" 
               "${PROJECT_SOURCE_DIR}/src/SiloReader.hh" 
         DESTINATION include/polytope)

# If we're parallel we have a few extra install items.
if (USE_MPI EQUAL 1)
  install (FILES "${PROJECT_SOURCE_DIR}/src/DistributedTessellator.hh" 
                 "${PROJECT_SOURCE_DIR}/src/SerialDistributedTessellator.hh" 
                 "${PROJECT_SOURCE_DIR}/src/checkDistributedTessellation.hh" 
           DESTINATION include/polytope)
endif ()

# build a CPack driven installer package
include(InstallRequiredSystemLibraries)
set (CPACK_RESOURCE_FILE_LICENSE  
     "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set (CPACK_PACKAGE_VERSION_MAJOR "${POLYTOPE_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${POLYTOPE_VERSION_MINOR}")
include(CPack)

