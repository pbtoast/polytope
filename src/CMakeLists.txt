# Include Silo stuff if we've found it.
if (SILO_FOUND)
  set(IO_SOURCES SiloWriter_2d.cc;SiloWriter_3d.cc)
endif()

# By default, we include the Voro tessellator(s).
set(TESSELLATOR_SOURCES VoroPP_2d.cc)
set(CWD ${CMAKE_CURRENT_SOURCE_DIR})

# Enable Triangle if triangle.h and triangle.c are found.
if (EXISTS ${CWD}/triangle.h AND 
    EXISTS ${CWD}/triangle.c)
  set(TESSELLATOR_SOURCES ${TESSELLATOR_SOURCES};TriangleTessellator.cc;triangle.c)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_TRIANGLE -DREAL=double -DTRILIBRARY -DANSI_DECLARATORS -DVOID=void")
endif ()

# Enable Tetgen if the proper files are found.
if (EXISTS ${CWD}/tetgen.h AND 
    EXISTS ${CWD}/tetgen.cxx AND 
    EXISTS ${CWD}/predicates.cxx)
  message("found tetgen!")
  set(TESSELLATOR_SOURCES ${TESSELLATOR_SOURCES};TetgenTessellator.cc;tetgen.cxx;predicates.cxx)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_TETGEN -DTETLIBRARY")
endif ()

add_library(polytope ${TESSELLATOR_SOURCES} ${IO_SOURCES})
target_link_libraries(polytope voro_2d)
