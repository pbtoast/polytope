# Serialization test.
add_executable(test_serialize test_serialize.cc)
add_test(test_serialize test_serialize)
target_link_libraries(test_serialize polytope)
set_tests_properties(test_serialize PROPERTIES PASS_REGULAR_EXPRESSION PASS
                     WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests)

# Convex hull 2D test.
add_executable(test_convexHull_2d test_convexHull_2d.cc)
add_test(test_convexHull_2d test_convexHull_2d)
target_link_libraries(test_convexHull_2d polytope)
set_tests_properties(test_convexHull_2d PROPERTIES PASS_REGULAR_EXPRESSION PASS
                     WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests)

# Convex hull 3D test.
add_executable(test_convexHull_3d test_convexHull_3d.cc)
add_test(test_convexHull_3d test_convexHull_3d)
target_link_libraries(test_convexHull_3d polytope)
set_tests_properties(test_convexHull_3d PROPERTIES PASS_REGULAR_EXPRESSION PASS
                     WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests)

# nearestPoint 2D test.
add_executable(test_nearestPoint_2d test_nearestPoint_2d.cc)
add_test(test_nearestPoint_2d test_nearestPoint_2d)
target_link_libraries(test_nearestPoint_2d polytope)
set_tests_properties(test_nearestPoint_2d PROPERTIES PASS_REGULAR_EXPRESSION PASS
                     WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests)

# Voro++ 2D tessellator.
add_executable(test1_VoroPP_2d test1_VoroPP_2d.cc)
add_test(test1_VoroPP_2d test1_VoroPP_2d)
target_link_libraries(test1_VoroPP_2d polytope)
set_tests_properties(test1_VoroPP_2d PROPERTIES PASS_REGULAR_EXPRESSION PASS
                     WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests)
if (Boost_FOUND)     
  add_executable(test_UnitSquare test_UnitSquare.cc)
  add_test(test_UnitSquare test_UnitSquare)
  target_link_libraries(test_UnitSquare polytope)
  set_tests_properties(test_UnitSquare PROPERTIES PASS_REGULAR_EXPRESSION PASS
                       WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests)
  add_executable(test_Degenerate test_Degenerate.cc)
  add_test(test_Degenerate test_Degenerate)
  target_link_libraries(test_Degenerate polytope)
  set_tests_properties(test_Degenerate PROPERTIES PASS_REGULAR_EXPRESSION PASS
                       WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests)
endif ()

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
  add_executable(test_RandomPoints test_RandomPoints.cc)
  add_test(test_RandomPoints test_RandomPoints)
  target_link_libraries(test_RandomPoints polytope ${SILO_LIBRARIES} ${HDF5_LIBRARIES})
  set_tests_properties(test_RandomPoints PROPERTIES PASS_REGULAR_EXPRESSION PASS
                       WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests)
  add_executable(test_PLCorder test_PLCorder.cc)
  add_test(test_PLCorder test_PLCorder)
  target_link_libraries(test_PLCorder polytope ${SILO_LIBRARIES} ${HDF5_LIBRARIES})
  set_tests_properties(test_PLCorder PROPERTIES PASS_REGULAR_EXPRESSION PASS
                       WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests)
  add_executable(test_Area test_Area.cc)
  add_test(test_Area test_Area)
  target_link_libraries(test_Area polytope ${SILO_LIBRARIES} ${HDF5_LIBRARIES})
  set_tests_properties(test_Area PROPERTIES PASS_REGULAR_EXPRESSION PASS
                       WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests)

  
  # Distributed Triangle tessellator.
  if (USE_MPI EQUAL 1)
    add_executable(test_DistributedTriangle test_DistributedTriangle.cc)
    add_test(test_DistributedTriangle test_DistributedTriangle)
    target_link_libraries(test_DistributedTriangle polytope)
    set_tests_properties(test_DistributedTriangle PROPERTIES PASS_REGULAR_EXPRESSION PASS
                         WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests)
    add_executable(test_SerialDistributedTriangle test_SerialDistributedTriangle.cc)
    add_test(test_SerialDistributedTriangle test_SerialDistributedTriangle)
    target_link_libraries(test_SerialDistributedTriangle polytope)
    set_tests_properties(test_SerialDistributedTriangle PROPERTIES PASS_REGULAR_EXPRESSION PASS
                         WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests)
    add_executable(test_FailedCommunication test_FailedCommunication.cc)
    add_test(test_FailedCommunication test_FailedCommunication)
    target_link_libraries(test_FailedCommunication polytope)
    set_tests_properties(test_FailedCommunication PROPERTIES PASS_REGULAR_EXPRESSION PASS
                         WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests)
    add_executable(test_DistributedRandomPoints test_DistributedRandomPoints.cc)
    add_test(test_DistributedRandomPoints test_DistributedRandomPoints)
    target_link_libraries(test_DistributedRandomPoints polytope)
    set_tests_properties(test_DistributedRandomPoints PROPERTIES PASS_REGULAR_EXPRESSION PASS
                         WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests)
    
  endif()
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

# Distributed Voro++ 3D tessellator.
  add_executable(test1_DistributedVoroPP_3d test1_DistributedVoroPP_3d.cc)
  add_test(test1_DistributedVoroPP_3d test1_DistributedVoroPP_3d)
  target_link_libraries(test1_DistributedVoroPP_3d polytope)
  set_tests_properties(test1_DistributedVoroPP_3d PROPERTIES PASS_REGULAR_EXPRESSION PASS
                       WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests)
endif()
