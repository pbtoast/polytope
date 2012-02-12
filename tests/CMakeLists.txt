# Convex hull 3D test.
add_executable(test_convexHull_3d test_convexHull_3d.cc)
add_test(test_convexHull_3d test_convexHull_3d)
target_link_libraries(test_convexHull_3d polytope)
set_tests_properties(test_convexHull_3d PROPERTIES PASS_REGULAR_EXPRESSION PASS
                     WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests)

# Voro++ 2D tessellator.
add_executable(test1_VoroPP_2d test1_VoroPP_2d.cc)
add_test(test1_VoroPP_2d test1_VoroPP_2d)
target_link_libraries(test1_VoroPP_2d polytope)
set_tests_properties(test1_VoroPP_2d PROPERTIES PASS_REGULAR_EXPRESSION PASS
                     WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests)

# Voro++ 3D tessellator.
add_executable(test1_VoroPP_3d test1_VoroPP_3d.cc)
add_test(test1_VoroPP_3d test1_VoroPP_3d)
target_link_libraries(test1_VoroPP_3d polytope)
set_tests_properties(test1_VoroPP_3d PROPERTIES PASS_REGULAR_EXPRESSION PASS
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

# Distributed Voro++ 2D tessellator.
if (USE_MPI EQUAL 1)
  add_executable(test1_DistributedVoroPP_2d test1_DistributedVoroPP_2d.cc)
  add_test(test1_DistributedVoroPP_2d test1_DistributedVoroPP_2d)
  target_link_libraries(test1_DistributedVoroPP_2d polytope)
  set_tests_properties(test1_DistributedVoroPP_2d PROPERTIES PASS_REGULAR_EXPRESSION PASS
                       WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests)
endif()
