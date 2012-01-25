add_executable(test1 test1.cc)

add_test(test1 test1)
set_tests_properties(test1 PROPERTIES PASS_REGULAR_EXPRESSION PASS
                     WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests)

