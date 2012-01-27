# Include Silo stuff if we've found it.
if (SILO_FOUND)
  set(IO_SOURCES SiloWriter_2d.cc;SiloWriter_3d.cc)
endif()

add_library(polytope VoroPP_2d.cc Triangle.cc ${IO_SOURCES})
target_link_libraries(polytope voro_2d)
