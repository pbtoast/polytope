# ----------------------------------------------------------------------------------
# CMake file for Polytope's Python bindings
#   -Auto-generates the pybindgen script PolytopeModulesBindings.py
#   -Auto-generates TessellatorModules.py based on the top-level config
#   -Builds and compiles the Python bindings according to UsePyBindGen cmake module
#
# NOTE: Modify the PYTHONPKGS variable to add/remove objects to the Python module
# ----------------------------------------------------------------------------------

# The objects we are adding to our Python module
set(PYTHONPKGS 
  "CXXContainers" 
  "MeshEditor" 
  "PLC" 
  "Tessellation" 
  "Tessellator" 
  "Functions"
)

# Generate the PolytopeModulesBindings.py driver script
configure_file(PolytopeModulesBindings.py.in  ${PYBINDGEN_DIR}/PolytopeModulesBindings.py)

# The macros for creating Python bindings
include(${PROJECT_SOURCE_DIR}/cmake/Modules/UsePyBindGen.cmake)

# The serial tessellators to be wrapped
set(SERIAL_TESSELLATORS "")
if(HAVE_TRIANGLE)
  set(SERIAL_TESSELLATORS ${SERIAL_TESSELLATORS} "TriangleTessellator2d")
endif()
if(HAVE_BOOST_VORONOI)
  set(SERIAL_TESSELLATORS ${SERIAL_TESSELLATORS} "BoostTessellator2d")
endif()
if(HAVE_TETGEN)
  set(SERIAL_TESSELLATORS ${SERIAL_TESSELLATORS} "TetgenTessellator3d")
endif()

# Generate the TessellatorModules input script
configure_file(TessellatorModule.py.in  ${PYBINDGEN_DIR}/TessellatorModule.py)

# Call the PyBindGen macros to create the C bindings
PYBINDGEN_GENERATE_BINDINGS("${PYTHONPKGS}")

# libpolytope location
set(POLYTOPE_STATIC_LIB "${PROJECT_BINARY_DIR}/src/libpolytope${CMAKE_STATIC_LIBRARY_SUFFIX}")

# libtriangle location
if(HAVE_TRIANGLE)
  set(TRIANGLE_STATIC_LIB "${PROJECT_BINARY_DIR}/src/libtriangle${CMAKE_STATIC_LIBRARY_SUFFIX}")
else()
  set(TRIANGLE_STATIC_LIB "")
endif()

# libtetgen location
if(HAVE_TETGEN)
  set(TETGEN_STATIC_LIB "${PROJECT_BINARY_DIR}/src/libtetgen${CMAKE_STATIC_LIBRARY_SUFFIX}")
else()
  set(TETGEN_STATIC_LIB "")
endif()

# libpython<version> location
string(REPLACE "bin/python" "lib" PYTHON_LINK_LIB ${PYTHON_EXE})
string(REPLACE "bin/python" "lib/libpython${PYTHON_VERSION}${CMAKE_STATIC_LIBRARY_SUFFIX}" PYTHON_STATIC_LIB ${PYTHON_EXE})
string(REPLACE "bin/python" "lib/libpython${PYTHON_VERSION}${CMAKE_SHARED_LIBRARY_SUFFIX}" PYTHON_SHARED_LIB ${PYTHON_EXE})

if (EXISTS "${PYTHON_STATIC_LIB}")
    set(PYTHON_LINK_LIB "${PYTHON_STATIC_LIB}")
else()
    set(PYTHON_LINK_LIB "${PYTHON_SHARED_LIB}")
endif()

# Include directories in the header file search path
include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_BINARY_DIR}/src")
include_directories("${PYTHON_INCLUDE_DIR}")
include_directories("${PYBINDGEN_DIR}")
include_directories("${CMAKE_INSTALL_PREFIX}/include")

set(PYTHON_LINK "python${PYTHON_VERSION}")

## Additional HDF5 link libraries
## NOTE:  This was added in to catch a few known build cases of static
##        silo and hdf5 on x86_64. Specifying the HDF5_ROOT will cause
##        Polytope to link against all the hdf5 libraries it finds.
## FIXME: A better solution is to allow users to specify a list of
##        additional link args when configuring Polytope.
#if (HDF5_ROOT)
#  set(HDF5_LINK "${HDF5_LIBRARIES}")
#else()
#  set(HDF5_LINK "")
#endif()

# On some Linux machines we get messed up with the HDF5 libs including 
# things like libdl and libm, so we try to cut down the list here to the 
# first couple of elements.
message("+++ ${HDF5_LIBRARIES}")
list(GET HDF5_LIBRARIES 0 1 HDF5_LIB)
#list(GET HDF5_LIBRARIES 0 3 HDF5_LIB)

link_directories("${PYTHON_LIB_DIR}")



# The static Python library of wrapped methods
add_library(PolytopeModulesStatic STATIC ${PYBINDGEN_GENERATED_SOURCE})
set_target_properties(PolytopeModulesStatic PROPERTIES OUTPUT_NAME PolytopeModules)

# The dynamic Python library of wrapped methods
add_library(PolytopeModules SHARED ${PYBINDGEN_GENERATED_SOURCE})
target_link_libraries(PolytopeModules 
#  ${PYTHON_LINK}
  ${POLYTOPE_STATIC_LIB} 
  ${TRIANGLE_STATIC_LIB}
  ${TETGEN_STATIC_LIB}
  ${SILO_LIBRARIES}
  ${HDF5_LIB}
)


#message("**** SYSTEM: ${CMAKE_SYSTEM}")
#message("**** CMAKE_SHARED_LINKER_FLAGS: ${CMAKE_SHARED_LINKER_FLAGS}")
#message("**** CMAKE_SHARED_LINKER_FLAGS_RELEASE: ${CMAKE_SHARED_LINKER_FLAGS_RELEASE}")
#message("**** CMAKE_CXX_IMPLICIT_LINK_LIBRARIES: ${CMAKE_CXX_IMPLICIT_LINK_LIBRARIES}")
#message("**** CMAKE_STANDARD_LIBRARIES: ${CMAKE_STANDARD_LIBRARIES}")
#message("**** CMAKE_LINK_DEPENDS_NO_SHARED: ${CMAKE_LINK_DEPENDS_NO_SHARED}")
#message("**** CMAKE_MODULE_LINKER_FLAGS: ${CMAKE_MODULE_LINKER_FLAGS}")
#message("**** CMAKE_MODULE_LINKER_FLAGS_RELEASE: ${CMAKE_MODULE_LINKER_FLAGS_RELEASE}")
#message("**** CMAKE_STANDARD_LIBRARIES: ${CMAKE_STANDARD_LIBRARIES}")
#message("**** HDF5_LIBRARIES: ${HDF5_LIBRARIES}")
#message("**** HDF5_C_LIBRARIES: ${HDF5_C_LIBRARIES}")
#message("**** HDF5_CXX_LIBRARIES: ${HDF5_CXX_LIBRARIES}")
#message("**** SILO_LIBRARIES: ${SILO_LIBRARIES}")
#set(CMAKE_VERBOSE_MAKEFILE TRUE)

# Remove the "lib" prefix convention from the library
# NOTE: "init<modulename>" is not defined properly if it begins with "lib"
set_target_properties(PolytopeModules PROPERTIES PREFIX "")

# Change .dylib suffix to .so on OSX
if (${APPLE})
  set_target_properties(PolytopeModules PROPERTIES SUFFIX ".so")

  # Force resolution of the python symbols to delay until we import the module.
  TARGET_LINK_LIBRARIES(PolytopeModules "-undefined suppress -flat_namespace")
endif()

# Allow rpath linkage to the Module
#set_target_properties(PolytopeModules PROPERTIES INSTALL_RPATH_USE_LINK_PATH true)

# Library install targets
install (TARGETS PolytopeModules DESTINATION lib/python${PYTHON_VERSION}/site-packages/polytope)
install (FILES polytope.pth DESTINATION lib/python${PYTHON_VERSION}/site-packages)
install (TARGETS PolytopeModulesStatic DESTINATION lib)
