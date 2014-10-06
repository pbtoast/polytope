#
# Change this list to add/remove modules wrapped by PyBindGen:
#
set(LIST_OF_WRAPPED_CLASSES "CXXContainers" "MeshEditor" "PLC" "Tessellation" "Tessellator" "Functions")

# The macros for creating Python bindings
include(${PROJECT_SOURCE_DIR}/cmake/Modules/UsePyBindGen.cmake)

# Additional args for Python execute
set(PYBINDGEN_ADDITIONAL_ARGS ${HAVE_TRIANGLE} ${HAVE_TETGEN} ${HAVE_BOOST_VORONOI} ${HAVE_MPI})

# Call the PyBindGen macros
PYBINDGEN_GENERATE_BINDINGS( "${LIST_OF_WRAPPED_CLASSES}" )

# libpolytope location
set(POLYTOPE_STATIC_LIB "${PROJECT_BINARY_DIR}/src/libpolytope${CMAKE_STATIC_LIBRARY_SUFFIX}")

# libtetgen location
if( TETGEN_FOUND )
  set(TETGEN_STATIC_LIB "${PROJECT_BINARY_DIR}/src/libtetgen.a")
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

set(PYTHON_LINK "python${PYTHON_VERSION}")

# On some Linux machines we get messed up with the HDF5 libs including things like libdl and libm, so 
# we try to cut down the list here to the first couple of elements.
list(GET HDF5_LIBRARIES 0 1 HDF5_LIB)

# The Python library of wrapped Polytope methods
link_directories("${PYTHON_LIB_DIR}")
add_library(PolytopeModules SHARED ${PYBINDGEN_GENERATED_SOURCE})
target_link_libraries(PolytopeModules 
#  ${PYTHON_LINK}
#  ${PYTHON_STATIC_LIB}
#  ${PYTHON_LINK_LIB}
  ${POLYTOPE_STATIC_LIB} 
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
