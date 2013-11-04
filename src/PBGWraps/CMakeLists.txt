# ----------------------------------------------------------------------------------
# CMake file for Polytope's Python bindings
#   -Auto-generates the pybindgen script PolytopeModulesBindings.py
#   -Auto-generates TessellatorModules.py based on the top-level config
#   -Builds and compiles the Python bindings according to UsePyBindGen cmake module
#
# NOTE: Modify the PYTHONPKGS variable to add/remove objects to the Python module
# ----------------------------------------------------------------------------------

# The objects we are adding to our Python module
set(PYTHONPKGS "CXXContainers" "MeshEditor" "PLC" "Tessellation" "Tessellator" "Functions")

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
set(POLYTOPE_STATIC_LIB "${PROJECT_BINARY_DIR}/src/libpolytope.a")

# libtetgen location
if(HAVE_TETGEN)
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


set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY TRUE)


# The static Python library of wrapped methods
add_library(PolytopeModules-static STATIC ${PYBINDGEN_GENERATED_SOURCE})
set_target_properties(PolytopeModules-static PROPERTIES OUTPUT_NAME PolytopeModules)
target_link_libraries(PolytopeModules-static 
  ${POLYTOPE_STATIC_LIB} ${TETGEN_STATIC_LIB} ${SILO_LIBRARIES})

# The dynamic Python library of wrapped methods
add_library(PolytopeModules SHARED ${PYBINDGEN_GENERATED_SOURCE})
target_link_libraries(PolytopeModules 
  ${POLYTOPE_STATIC_LIB} ${TETGEN_STATIC_LIB} ${SILO_LIBRARIES})

# Remove the "lib" prefix convention from the library
# NOTE: "init<modulename>" is not defined properly if it begins with "lib"
set_target_properties(PolytopeModules PROPERTIES PREFIX "")

# Allow rpath linkage to the Module
set_target_properties(PolytopeModules PROPERTIES INSTALL_RPATH_USE_LINK_PATH false)


# Library install targets
install (TARGETS PolytopeModules DESTINATION lib)
#install (TARGETS PolytopeModules DESTINATION lib OPTIONAL)
#install (TARGETS PolytopeModules-static DESTINATION lib OPTIONAL)
