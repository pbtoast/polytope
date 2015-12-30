# Include Silo stuff if we've found it.
if (SILO_FOUND)
  set(IO_SOURCES SiloWriter_2d.cc SiloWriter_3d.cc SiloReader_2d.cc SiloReader_3d.cc SiloUtils.cc)
endif()

# By default, we include the Voro tessellator(s).
set(TESSELLATOR_SOURCES VoroPP_2d.cc VoroPP_3d.cc KeyTraits.cc predicates.cc PLC_CSG.cc)
set(CWD ${CMAKE_CURRENT_SOURCE_DIR})

# Check if we can include each of the serial tessellators:
if (USE_BOOST)
  list(APPEND TESSELLATOR_SOURCES BoostOrphanage.cc MeshEditor.cc)
endif()

# TriangleTessellator and triangle source
if (HAVE_TRIANGLE)
  list(APPEND TESSELLATOR_SOURCES TriangleTessellator.cc)
  add_library(triangle STATIC  predicates.cc triangle.c)
  #set(TRIANGLE_LIB triangle PARENT_SCOPE)
  set(TRIANGLE_LIB triangle)
  set_target_properties(triangle PROPERTIES
    COMPILE_FLAGS "${COMPILE_FLAGS} -DTRILIBRARY -DCDT_ONLY -DANSI_DECLARATIONS -fPIC -DNO_TIMER"
    LINKER_LANGUAGE C)
endif()

# BoostTessellator
if (HAVE_BOOST_VORONOI)
  list(APPEND TESSELLATOR_SOURCES BoostTessellator.cc)
endif ()

# TetgenTessellator and tetgen source
if (HAVE_TETGEN)
  list(APPEND TESSELLATOR_SOURCES TetgenTessellator.cc)
  add_library(tetgen tetgen.cxx predicates.cc)
  #set(TETGEN_LIB tetgen PARENT_SCOPE)
  set(TETGEN_LIB tetgen)
  set_target_properties(${TETGEN_LIB} PROPERTIES 
    COMPILE_FLAGS "${COMPILE_FLAGS} -DTETLIBRARY -Wuninitialized")
endif ()

# Check if we can incude the distributed tessellators:
if (HAVE_MPI)
  list(APPEND TESSELLATOR_SOURCES DistributedTessellator.cc SerialDistributedTessellator.cc)
endif()

# Add header files for IDE projects
file(GLOB HEADERS *.hh)
list(APPEND HEADERS "${PROJECT_BINARY_DIR}/polytope.hh")

# On some Linux machines we get messed up with the HDF5 libs including 
# things like libdl and libm, so we try to cut down the list here to the 
# first couple of elements.
list(GET HDF5_LIBRARIES 0 1 HDF5_LIB)

add_library(polytope STATIC
            ${TESSELLATOR_SOURCES}
	    ${IO_SOURCES}
	    ErrorHandler.cc
	    polytope_internal_abort.cc
	    ${HEADERS})

target_link_libraries(polytope
                      voro_2d 
		      voro_3d 
		      ${TRIANGLE_LIB} 
		      ${TETGEN_LIB} 
		      ${SILO_LIBRARIES}
		      ${HDF5_LIBRARIES}
		      ${MPI_C_LIBRARIES})

# Library install targets
install(TARGETS polytope ${TRIANGLE_LIB} ${TETGEN_LIB} DESTINATION lib)
