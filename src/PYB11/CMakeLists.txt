# ----------------------------------------------------------------------------------
# CMake file for Polytope's Python bindings using PYB11Generator
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
configure_file(PolytopeModulesBindings.py.in ${PYBINDGEN_DIR}/PolytopeModulesBindings.py)

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
#configure_file(TessellatorModule.py.in  ${PYBINDGEN_DIR}/TessellatorModule.py)
configure_file(TessellatorModule.py.in  ${PROJECT_BINARY_DIR}/src/PBGWraps/TessellatorModule.py)

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

# On some Linux machines we get messed up with the HDF5 libs including 
# things like libdl and libm, so we try to cut down the list here to the 
# first couple of elements.
list(GET HDF5_LIBRARIES 0 1 HDF5_LIB)

link_directories("${PYTHON_LIB_DIR}")

# The static Python library of wrapped methods
add_library(PolytopeModulesStatic STATIC ${PYBINDGEN_GENERATED_SOURCE})
set_target_properties(PolytopeModulesStatic PROPERTIES OUTPUT_NAME PolytopeModules)

# The dynamic Python library of wrapped methods
add_library(PolytopeModules SHARED ${PYBINDGEN_GENERATED_SOURCE})
target_link_libraries(PolytopeModules 
                      #${PYTHON_LINK}
                      ${POLYTOPE_STATIC_LIB} 
		      ${TRIANGLE_STATIC_LIB}
		      ${TETGEN_STATIC_LIB})

if (HAVE_SILO)
  target_link_libraries(PolytopeModules 
                        ${SILO_LIBRARIES})
endif()
if (HAVE_HDF5)
  target_link_libraries(PolytopeModules 
                        ${HDF5_LIB})
endif()

# Remove the "lib" prefix convention from the library
# NOTE: "init<modulename>" is not defined properly if it begins with "lib"
set_target_properties(PolytopeModules PROPERTIES PREFIX "")

# Change .dylib suffix to .so on OSX
if (${APPLE})
  set_target_properties(PolytopeModules PROPERTIES SUFFIX ".so")

  # Force resolution of the python symbols to delay until we import the module.
  target_link_libraries(PolytopeModules "-undefined suppress -flat_namespace")
endif()

# Library install targets
# We can't count on installing in the PYTHON_SITE_PACKAGES when python is a system place.
#if (PYTHON_EXE)
#  install(FILES polytope.pth DESTINATION ${PYTHON_SITE_PACKAGES})
#  install(TARGETS PolytopeModules DESTINATION ${PYTHON_SITE_PACKAGES}/polytope)
#else()
  install(FILES polytope.pth DESTINATION lib/python${PYTHON_VERSION}/site-packages)
  install(TARGETS PolytopeModules DESTINATION lib/python${PYTHON_VERSION}/site-packages/polytope)
#endif()

install(TARGETS PolytopeModulesStatic DESTINATION lib)



