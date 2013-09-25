add_library(polytope_c ErrorHandler.cc polytope_internal_abort.cc)
target_link_libraries(polytope_c polytope)

# Library install targets
install (TARGETS polytope_c DESTINATION lib)
