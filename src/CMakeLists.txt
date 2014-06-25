# Include Silo stuff if we've found it.
if (SILO_FOUND)
  set(IO_SOURCES SiloWriter_2d.cc SiloWriter_3d.cc SiloReader_2d.cc SiloReader_3d.cc SiloUtils.cc)
  message("++ IO_SOURCES: ${IO_SOURCES}")
  message("++ SILO_LIBRARIES: ${SILO_LIBRARIES}")
  message("++ HDF5_LIBRARIES: ${HDF5_LIBRARIES}")    
endif()

# By default, we include the Voro tessellator(s).
set(TESSELLATOR_SOURCES VoroPP_2d.cc VoroPP_3d.cc KeyTraits.cc predicates.cc PLC_CSG.cc)
set(CWD ${CMAKE_CURRENT_SOURCE_DIR})

if (HAVE_BOOST)
  set(TESSELLATOR_SOURCES ${TESSELLATOR_SOURCES} BoostOrphanage.cc MeshEditor.cc)
endif()

if (TRIANGLE_FOUND)
  set(TESSELLATOR_SOURCES ${TESSELLATOR_SOURCES} TriangleTessellator.cc)
endif()

if (HAVE_BOOST_VORONOI)
  set(TESSELLATOR_SOURCES ${TESSELLATOR_SOURCES} BoostTessellator.cc)
endif ()

if (TETGEN_FOUND)
  set(TESSELLATOR_SOURCES ${TESSELLATOR_SOURCES} TetgenTessellator.cc)
  add_library(tetgen tetgen.cxx predicates.cc)
  set(TETGEN_LIB tetgen)
  set_target_properties(${TETGEN_LIB} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DTETLIBRARY")
else ()
  set(TETGEN_LIB "")
endif ()

if (USE_MPI EQUAL 1)
  set(TESSELLATOR_SOURCES ${TESSELLATOR_SOURCES} DistributedTessellator.cc 
                                                 SerialDistributedTessellator.cc)
endif()

add_library(polytope ${TESSELLATOR_SOURCES} ${IO_SOURCES} ErrorHandler.cc polytope_internal_abort.cc)
target_link_libraries(polytope voro_2d voro_3d ${TETGEN_LIB} ${SILO_LIBRARIES} ${HDF5_LIBRARIES})

# Library install targets
install (TARGETS polytope ${TETGEN_LIB} DESTINATION lib)
