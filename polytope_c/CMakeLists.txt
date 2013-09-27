# Build a configuration header file from our options.

configure_file(
  "${PROJECT_SOURCE_DIR}/polytope_c/polytope_c.h.in"
  "${PROJECT_BINARY_DIR}/polytope_c.h"
)

add_library(polytope_c polytope_plc.cc polytope_tessellator.cc 
            polytope_read_silo.cc polytope_write_silo.cc)
target_link_libraries(polytope_c polytope)

# Library install targets
install (TARGETS polytope_c DESTINATION lib)
