# ------------------------------------------------------------------#
# Change this list to add/remove modules wrapped by PyBindGen:      #
# ------------------------------------------------------------------#
set(LIST_OF_WRAPPED_CLASSES 
  "CXXContainers" 
  "MeshEditor" 
  "PLC" 
  "Tessellation" 
  "Tessellator" 
  "Functions"
)


# The macros for creating Python bindings
include(${PROJECT_SOURCE_DIR}/cmake/Modules/UsePyBindGen.cmake)

# Additional args for Python execute
set(PYBINDGEN_ADDITIONAL_ARGS ${HAVE_TRIANGLE} ${HAVE_TETGEN} ${HAVE_BOOST_VORONOI} ${HAVE_MPI})

# Call the PyBindGen macros
PYBINDGEN_GENERATE_BINDINGS( "${LIST_OF_WRAPPED_CLASSES}" )

# libpolytope location
set(POLYTOPE_STATIC_LIB "${PROJECT_BINARY_DIR}/src/libpolytope.a")

# libtetgen location
if( TETGEN_FOUND )
  set(TETGEN_STATIC_LIB "${PROJECT_BINARY_DIR}/src/libtetgen.a")
endif()

# libpython<version> location
string(REPLACE "bin/python" "lib" PYTHON_LINK_LIB ${PYTHON_EXE})
string(REPLACE "bin/python" "lib/libpython${PYTHON_VERSION}.a"
  PYTHON_STATIC_LIB ${PYTHON_EXE})

# Include directories in the header file search path
include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_BINARY_DIR}/src")
include_directories("${PYTHON_INCLUDE_DIR}")
include_directories("${PYBINDGEN_DIR}")

set(PYTHON_LINK "python")

# The Python library of wrapped Polytope methods
add_library(PolytopeModules SHARED ${PYBINDGEN_GENERATED_SOURCE})

# Change .dylib suffix to .so on OSX and add additional python linking
if (${APPLE})
  set(PYTHON_LINK "python")
  set_target_properties(PolytopeModules PROPERTIES SUFFIX ".so")
else()
  set(PYTHON_LINK "")
endif()

# Additional HDF5 link libraries
# NOTE:  This was added in to catch a few known build cases of static
#        silo and hdf5 on x86_64. Specifying the HDF5_ROOT will cause
#        Polytope to link against all the hdf5 libraries it finds.
# FIXME: A better solution is to allow users to specify a list of
#        additional link args when configuring Polytope.
if (HDF5_ROOT)
  set(HDF5_LINK "${HDF5_LIBRARIES}")
else()
  set(HDF5_LINK "")
endif()

# Remove the "lib" prefix convention from the library so Python can import
set_target_properties(PolytopeModules PROPERTIES PREFIX "")

# Allow rpath linkage to the Module
# NOTE: Commented out rpath linking during portability testing.
#       Haven't tested if turning back on will cause problems.
#set_target_properties(PolytopeModules PROPERTIES INSTALL_RPATH_USE_LINK_PATH true)

target_link_libraries(PolytopeModules 
  ${PYTHON_LINK}
#  ${PYTHON_STATIC_LIB}
  ${POLYTOPE_STATIC_LIB} 
  ${TETGEN_STATIC_LIB}
  ${SILO_LIBRARIES}
  ${HDF5_LINK}
)

# Remove the "lib" prefix convention from the library
# NOTE: "init<modulename>" is not defined properly if it begins with "lib"
set_target_properties(PolytopeModules PROPERTIES PREFIX "")

# Change .dylib suffix to .so on OSX
if (${APPLE})
  set_target_properties(PolytopeModules PROPERTIES SUFFIX ".so")
endif()

# Allow rpath linkage to the Module
#set_target_properties(PolytopeModules PROPERTIES INSTALL_RPATH_USE_LINK_PATH true)

# Library install targets
install (TARGETS PolytopeModules DESTINATION lib)
