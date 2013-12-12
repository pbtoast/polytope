# Include Silo stuff if we've found it.
if (SILO_FOUND)
  set(IO_SOURCES SiloWriter_2d.cc SiloWriter_3d.cc SiloReader_2d.cc SiloReader_3d.cc SiloUtils.cc)
endif()

# By default, we include the Voro tessellator(s).
set(TESSELLATOR_SOURCES VoroPP_2d.cc VoroPP_3d.cc KeyTraits.cc predicates.cc PLC_CSG.cc MeshEditor.cc)
set(CWD ${CMAKE_CURRENT_SOURCE_DIR})

# Check if we can include each of the serial tessellators:
if (HAVE_BOOST)
  set(TESSELLATOR_SOURCES ${TESSELLATOR_SOURCES} BoostOrphanage.cc)
endif()

if (HAVE_TRIANGLE)
  set(TESSELLATOR_SOURCES ${TESSELLATOR_SOURCES} TriangleTessellator.cc)
endif()

if (HAVE_BOOST_VORONOI)
  set(TESSELLATOR_SOURCES ${TESSELLATOR_SOURCES} BoostTessellator.cc)
endif ()

if (HAVE_TETGEN)
  set(TESSELLATOR_SOURCES ${TESSELLATOR_SOURCES} TetgenTessellator.cc)
  add_library(tetgen tetgen.cxx predicates.cc)
  set(TETGEN_LIB tetgen)
  set_target_properties(${TETGEN_LIB} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DTETLIBRARY")
else ()
  set(TETGEN_LIB "")
endif ()

# Check if we can incude the distributed tessellators:
if (HAVE_MPI)
  set(TESSELLATOR_SOURCES ${TESSELLATOR_SOURCES} DistributedTessellator.cc 
                                                 SerialDistributedTessellator.cc)
endif()

# The Polytope library:
add_library(polytope ${TESSELLATOR_SOURCES} ${IO_SOURCES} ErrorHandler.cc polytope_internal_abort.cc)
target_link_libraries(polytope voro_2d voro_3d ${TETGEN_LIB} ${SILO_LIBRARIES} ${HDF5_LIBRARIES})

# Library install targets
install (TARGETS polytope ${TETGEN_LIB} DESTINATION lib)
