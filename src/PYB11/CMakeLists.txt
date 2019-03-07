# ----------------------------------------------------------------------------------
# CMake file for Polytope's Python bindings using PYB11Generator
#   -Builds and compiles the Python bindings according to UsePyBindGen cmake module
# ----------------------------------------------------------------------------------

# The PYB11Generator module file.
set(PYB11_MODULE_NAME "polytope")

#
set (PYB11_SRC
     PLC.py
     ReducedPLC.py
     Tessellation.py)

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

# Generate the Tessellators.py script with the tessellators we're building
configure_file(Tessellators.py.in ${PROJECT_BINARY_DIR}/src/PYB11/Tessellators.py)

# The macros for creating Python bindings
include(${PROJECT_SOURCE_DIR}/cmake/Modules/UsePYB11Generator.cmake)

# Use PYB11Generator to create the pybind11 bindings
PYB11_GENERATE_BINDINGS()

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
include_directories("${PYB11_DIR}")
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
add_library(pypolytopeStatic STATIC ${PYB11_GENERATED_SOURCE})
set_target_properties(pypolytopeStatic PROPERTIES OUTPUT_NAME polytope)

# The dynamic Python library of wrapped methods
add_library(pypolytope SHARED ${PYB11_GENERATED_SOURCE})
set_target_properties(pypolytope PROPERTIES OUTPUT_NAME polytope)
target_link_libraries(pypolytope 
                      #${PYTHON_LINK}
                      ${POLYTOPE_STATIC_LIB} 
		      ${TRIANGLE_STATIC_LIB}
		      ${TETGEN_STATIC_LIB})

if (HAVE_SILO)
  target_link_libraries(pypolytope 
                        ${SILO_LIBRARIES})
endif()
if (HAVE_HDF5)
  target_link_libraries(pypolytope 
                        ${HDF5_LIB})
endif()

# Remove the "lib" prefix convention from the library
# NOTE: "init<modulename>" is not defined properly if it begins with "lib"
set_target_properties(pypolytope PROPERTIES PREFIX "")

# Change .dylib suffix to .so on OSX
if (${APPLE})
  set_target_properties(pypolytope PROPERTIES SUFFIX ".so")

  # Force resolution of the python symbols to delay until we import the module.
  target_link_libraries(pypolytope "-undefined suppress -flat_namespace")
endif()

# Library install targets
# We can't count on installing in the PYTHON_SITE_PACKAGES when python is a system place.
install(FILES polytope.pth DESTINATION lib/python${PYTHON_VERSION}/site-packages)
install(TARGETS pypolytope DESTINATION lib/python${PYTHON_VERSION}/site-packages/polytope)

install(TARGETS pypolytopeStatic DESTINATION lib)
