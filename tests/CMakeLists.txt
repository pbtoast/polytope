add_executable(test1 test1.cc)
add_executable(test1_VoroPP_2d test1_VoroPP_2d.cc)

add_test(test1 test1)
set_tests_properties(test1 PROPERTIES PASS_REGULAR_EXPRESSION PASS
                     WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests)

add_test(test1_VoroPP_2d test1_VoroPP_2d)
target_link_libraries(test1_VoroPP_2d polytope)
set_tests_properties(test1_VoroPP_2d PROPERTIES PASS_REGULAR_EXPRESSION PASS
                     WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests)
