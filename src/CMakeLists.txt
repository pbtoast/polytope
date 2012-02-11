# Include Silo stuff if we've found it.
if (SILO_FOUND)
  set(IO_SOURCES SiloWriter_2d.cc SiloWriter_3d.cc)
endif()

# By default, we include the Voro tessellator(s).
set(TESSELLATOR_SOURCES VoroPP_2d.cc VoroPP_3d.cc predicates.cc)
set(CWD ${CMAKE_CURRENT_SOURCE_DIR})

if (TRIANGLE_FOUND)
  set(TESSELLATOR_SOURCES ${TESSELLATOR_SOURCES} TriangleTessellator.cc)
endif()

if (TETGEN_FOUND)
  set(TESSELLATOR_SOURCES ${TESSELLATOR_SOURCES} TetgenTessellator.cc)
endif ()

if (USE_MPI EQUAL 1)
  set(TESSELLATOR_SOURCES ${TESSELLATOR_SOURCES} DistributedTessellator.cc)
endif()

add_library(polytope ${TESSELLATOR_SOURCES} ${IO_SOURCES} ErrorHandler.cc)
target_link_libraries(polytope voro_2d voro_3d ${SILO_LIBRARIES} ${HDF5_LIBRARIES})

# Library install targets
install (TARGETS polytope DESTINATION lib)
