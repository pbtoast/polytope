# Voro++ 2D tessellator.
add_executable(test1_VoroPP_2d test1_VoroPP_2d.cc)
add_test(test1_VoroPP_2d test1_VoroPP_2d test1_TriangleTessellator)
target_link_libraries(test1_VoroPP_2d polytope)
set_tests_properties(test1_VoroPP_2d PROPERTIES PASS_REGULAR_EXPRESSION PASS
                     WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests)

# Triangle (2D) tessellator.
if (TRIANGLE_FOUND)
  add_executable(test_TriangleTessellator test_TriangleTessellator.cc)
  add_test(test_TriangleTessellator test_TriangleTessellator)
  target_link_libraries(test_TriangleTessellator polytope ${SILO_LIBRARIES} ${HDF5_LIBRARIES})
  set_tests_properties(test_TriangleTessellator PROPERTIES PASS_REGULAR_EXPRESSION PASS
                       WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests)
endif()

# Tetgen (3D) tessellator.
if (TETGEN_FOUND)
  add_executable(test_TetgenTessellator test_TetgenTessellator.cc)
  add_test(test_TetgenTessellator test_TetgenTessellator)
  target_link_libraries(test_TetgenTessellator polytope ${SILO_LIBRARIES} ${HDF5_LIBRARIES})
  set_tests_properties(test_TetgenTessellator PROPERTIES PASS_REGULAR_EXPRESSION PASS
                       WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests)
endif()
