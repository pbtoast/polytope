# Minimum CMake version.
cmake_minimum_required (VERSION 2.8.5)

project(Polytope)

##set(CMAKE_SKIP_RPATH ON)

# Project version numbers.
set(POLYTOPE_VERSION_MAJOR 1)
set(POLYTOPE_VERSION_MINOR 0)

# Adjust Polytope's module path.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Set compilers. This must be done before enabling languages.
if (CC)
  set(CMAKE_C_COMPILER "${CC}")
  message(STATUS "C compiler is ${CMAKE_C_COMPILER}")
endif()
if (CXX)
  set(CMAKE_CXX_COMPILER "${CXX}")
  message(STATUS "C++ compiler is ${CMAKE_CXX_COMPILER}")
endif()

# Build everything as static libs.
option(BUILD_SHARED_LIBS "Build Polytope as shared library" OFF)

# Figure out the system type.
if (APPLE)
  set(SYS_FLAGS "-DAPPLE=1")
else ()
  if (LINUX)
  set(SYS_FLAGS "-DLINUX=1")
  endif ()
endif ()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SYS_FLAGS}")

# General compiler flags.
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-sign-compare -ansi -fPIC")
endif()
if (CMAKE_CXX_COMPILER_ID STREQUAL "XL")
  set(USING_XLC 1)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qPIC")
else()
  if (NOT MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  endif()
endif()
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-sign-compare")
endif()
if (MSVC)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS -D_USE_MATH_DEFINES)
endif (MSVC)

# Figure out MPI.
option(USE_MPI "Compile MPI parallel tessellator" OFF)
if (USE_MPI)
  # CC and CXX should already have been set in Makefile or wherever.
  find_package(MPI REQUIRED)
  set(HAVE_MPI 1)
  include_directories("${MPI_CXX_INCLUDE_PATH}")
endif()

# Find Boost.
option(USE_BOOST "Use Boost Voronoi tessellator" ON)

### Disable this on IBM xlC since Boost.Geometry does not compile correctly there.
##if (USE_BOOST AND USING_XLC)
##  message("-- Boost being disabled on IBM for now due to incompatibility with xlC")
##  set(USE_BOOST OFF)0
##endif ()

if (USE_BOOST)
  set(Boost_ADDITIONAL_VERSIONS "1.51" "1.51.0" "1.52" "1.52.0" "1.53" "1.53.0" "1.60")
  find_package(Boost 1.50 REQUIRED)
  if (Boost_FOUND)
    include_directories("${Boost_INCLUDE_DIR}")
    set(HAVE_BOOST 1)
    message(STATUS "Found Boost installation at ${Boost_INCLUDE_DIR}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBoost_FOUND=1")
    if (Boost_MINOR_VERSION GREATER 51 AND
	EXISTS ${Boost_INCLUDE_DIR}/boost/polygon/voronoi.hpp)
      set(HAVE_BOOST_VORONOI 1)
      message(STATUS "Found Boost Voronoi in ${Boost_INCLUDE_DIR}/boost/polygon")
    else ()
      message(STATUS "Boost Voronoi not found. BoostTessellator is disabled.")
    endif ()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBoost_FOUND=1")
  endif ()
endif ()

# Find Python
option(USE_PYTHON "Generate Python interace" OFF)
if (USE_PYTHON)
  
  if (PYTHON_EXE)
    set(PYTHON_EXECUTABLE ${PYTHON_EXE})
    include(${PROJECT_SOURCE_DIR}/cmake/Modules/FindPythonVersion.cmake)
    find_python_version()
    string(REPLACE "/bin/python" "" PYTHON_ROOT ${PYTHON_EXE})
    set(PYTHON_INCLUDE_DIR "${PYTHON_ROOT}/include/python${PYTHON_VERSION}")
    set(PYTHON_LIB_DIR "${PYTHON_ROOT}/lib")
    set(PYTHON_LIBRARY "${PYTHON_LIB_DIR}/libpython.a")
  else()
    find_package(PythonLibs REQUIRED)
    find_package(PythonInterp REQUIRED)
    set(PYTHON_VERSION "${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}")
  endif()
  
  if (PYTHON_EXECUTABLE)
    include(${PROJECT_SOURCE_DIR}/cmake/Modules/FindPythonModule.cmake)
    find_python_module(pybindgen REQUIRED)
    set(PYBINDGEN_DIR "${PROJECT_SOURCE_DIR}/src/PBGWraps")
    set(PYBINDGEN_MODULE_NAME "PolytopeModules")
    if (PY_PYBINDGEN)
      set(HAVE_PYTHON 1)
      message(STATUS "Found Python version ${PYTHON_VERSION}")
      message(STATUS "Using Python executable ${PYTHON_EXECUTABLE}")
    endif()
  endif()

  if (HAVE_PYTHON)
    include(${PROJECT_SOURCE_DIR}/cmake/Modules/FindPythonSitePackages.cmake)
    find_python_site_packages()
  endif()
endif()

#Check for HDF5 and Silo
option(USE_SILO "Use SILO reader/writer" OFF)
if (USE_SILO)
  # Find HDF5.
  set(HDF5_USE_STATIC_LIBRARIES 1)
  set(ENV{HDF5_ROOT} ${HDF5_ROOT})
  find_package(HDF5 REQUIRED)
  if (HDF5_FOUND)
    set(HAVE_HDF5 1)
    # Get rid of'debug' and 'optimized' detritis.
    list(REMOVE_ITEM HDF5_LIBRARIES debug optimized)
    list(REMOVE_ITEM HDF5_INCLUDE_DIRS debug optimized)

    # Get rid of duplicates
    # NOTE: Taken from Will Dicharry's FindHDF5.cmake module. This ensures
    #       duplicates are removed from the beginning of the list so as
    #       not to mess with the unresolved symbols search at linking.
    list(REVERSE HDF5_LIBRARIES)
    list(REMOVE_DUPLICATES HDF5_LIBRARIES)
    list(REVERSE HDF5_LIBRARIES)
    message(STATUS "Found HDF5 files in ${HDF5_INCLUDE_DIRS}")

    # Do we have parallel HDF5?
    if (USE_MPI)
      if (HDF5_IS_PARALLEL)
        set(HAVE_PARALLEL_HDF5 1)
      endif ()
    endif ()
  endif()

  find_package(Silo REQUIRED)
  if (SILO_FOUND)
    set(HAVE_SILO 1)
    include_directories("${SILO_INCLUDE_DIRS}")
    include_directories("${HDF5_INCLUDE_DIRS}")
#         FIXME: For now, we assume libsiloh5 lives with libhdf5.
#         link_directories("${SILO_LIBRARY_DIRS}")
    link_directories("${HDF5_LIBRARY_DIRS}")
  endif()
endif()

if(HDF5_LIBRARIES)
  _remove_duplicates_from_beginning(HDF5_LIBRARIES)
endif()


# Enable Triangle if triangle.h and triangle.c are found.
if (EXISTS ${PROJECT_SOURCE_DIR}/src/triangle.h AND 
    EXISTS ${PROJECT_SOURCE_DIR}/src/triangle.c AND
    Boost_FOUND)
  set(HAVE_TRIANGLE 1)
  message(STATUS "Triangle has been found! TriangleTessellator is enabled.")
else()
  message(STATUS "Triangle not found! TriangleTessellator is disabled.")
endif ()

# Enable Tetgen if the proper files are found.
if (EXISTS ${PROJECT_SOURCE_DIR}/src/tetgen.h AND 
    EXISTS ${PROJECT_SOURCE_DIR}/src/tetgen.cxx)  
  set(HAVE_TETGEN 1)
  message(STATUS "Tetgen has been found! TetgenTessellator is enabled.")
else()
  message(STATUS "Tetgen not found! TetgenTessellator is disabled.")
endif ()

# Build a configuration header file from our options.
configure_file(
  "${PROJECT_SOURCE_DIR}/src/polytope.hh.in"
  "${PROJECT_BINARY_DIR}/polytope.hh"
)

option(HEADER_ONLY "Header only installation" OFF)
if (NOT HEADER_ONLY)
  # Include the binary directory in the header file search path.
  include_directories("${PROJECT_BINARY_DIR}")
  include_directories("${PROJECT_SOURCE_DIR}/src")

  # Libraries

  # # Voro 2d/3d libraries.
  #include_directories(voro_2d)
  #add_subdirectory(voro_2d)
  #include_directories(voro_3d)
  #add_subdirectory(voro_3d)

  # Polytope proper.
  include_directories("${PROJECT_SOURCE_DIR}/src/")
  add_subdirectory(src)

  # Tests directory
  option(TESTING "Generate tests" ON)
  if (TESTING)
    enable_testing()
    add_subdirectory(tests)
    include(CTest)
  endif ()
  
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
	       "${PROJECT_SOURCE_DIR}/src/MeshEditor.hh"
               "${PROJECT_SOURCE_DIR}/src/ErrorHandler.hh" 
               "${PROJECT_SOURCE_DIR}/src/PLC.hh" 
	       "${PROJECT_SOURCE_DIR}/src/PLC_CSG_2d.hh"
	       "${PROJECT_SOURCE_DIR}/src/PLC_CSG_3d.hh"
               "${PROJECT_SOURCE_DIR}/src/ReducedPLC.hh"
               "${PROJECT_SOURCE_DIR}/src/simplifyPLCfacets.hh"
               "${PROJECT_SOURCE_DIR}/src/polytope_internal.hh" 
               "${PROJECT_SOURCE_DIR}/src/QuantTessellation.hh" 
	       "${PROJECT_SOURCE_DIR}/src/QuantizedCoordinates.hh"
               "${PROJECT_SOURCE_DIR}/src/VoroPP_2d.hh" 
               "${PROJECT_SOURCE_DIR}/src/VoroPP_3d.hh" 
               "${PROJECT_SOURCE_DIR}/src/TriangleTessellator.hh"
	       "${PROJECT_SOURCE_DIR}/src/BoostTessellator.hh"
	       "${PROJECT_SOURCE_DIR}/src/BoostTessellatorTraits.hh"
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
               "${PROJECT_SOURCE_DIR}/src/polytope_serialize.hh" 
               "${PROJECT_SOURCE_DIR}/src/polytope_geometric_utilities.hh" 
               "${PROJECT_SOURCE_DIR}/src/polytope_parallel_utilities.hh"
	       "${PROJECT_SOURCE_DIR}/src/polytope_tessellator_utilities.hh"
	       "${PROJECT_SOURCE_DIR}/src/PLC_CSG_2d.hh"
	       "${PROJECT_SOURCE_DIR}/src/PLC_CSG_3d.hh"
	       "${PROJECT_SOURCE_DIR}/src/PLC_Boost_2d.hh"
               "${PROJECT_SOURCE_DIR}/src/SiloWriter.hh" 
               "${PROJECT_SOURCE_DIR}/src/SiloReader.hh" 
               "${PROJECT_SOURCE_DIR}/src/polytope_write_OOGL.hh" 
         DESTINATION include/polytope)

# If we're parallel we have a few extra install items.
if (USE_MPI)
  install(FILES "${PROJECT_SOURCE_DIR}/src/DistributedTessellator.hh" 
                "${PROJECT_SOURCE_DIR}/src/SerialDistributedTessellator.hh" 
                "${PROJECT_SOURCE_DIR}/src/checkDistributedTessellation.hh" 
          DESTINATION include/polytope)
endif ()

# build a CPack driven installer package
include(InstallRequiredSystemLibraries)
set(CPACK_RESOURCE_FILE_LICENSE  
  "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_VERSION_MAJOR "${POLYTOPE_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${POLYTOPE_VERSION_MINOR}")
include(CPack)
