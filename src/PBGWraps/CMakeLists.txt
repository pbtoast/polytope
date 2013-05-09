#
# Change this list to add/remove modules wrapped by PyBindGen:
#
set(LIST_OF_WRAPPED_CLASSES "CXXContainers" "MeshEditor" "PLC" "Tessellation" "Tessellator")

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

# The Python library of wrapped Polytope methods
add_library(PolytopeModules SHARED ${PYBINDGEN_GENERATED_SOURCE})
target_link_libraries(PolytopeModules 
  ${POLYTOPE_STATIC_LIB} 
  ${TETGEN_STATIC_LIB}
#  ${PYTHON_STATIC_LIB}
)

# Remove the "lib" prefix convention from the library
# NOTE: "init<modulename>" is not defined properly if it begins with "lib"
set_target_properties(PolytopeModules PROPERTIES PREFIX "")

# Allow rpath linkage to the Module
#set_target_properties(PolytopeModules PROPERTIES INSTALL_RPATH_USE_LINK_PATH true)

# Library install targets
install (TARGETS PolytopeModules DESTINATION lib)
