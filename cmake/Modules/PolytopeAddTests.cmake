# CMake macro for adding tests into Polytope
# Consult test/CMakeLists.txt for usage

macro(POLYTOPE_ADD_TEST name dependency_list)
  set(BUILD_TEST true)
  foreach(_dependency ${dependency_list})
    set(DEP_NAME "HAVE_${_dependency}")
    if(NOT ${DEP_NAME})
      set(BUILD_TEST false)
    endif()
  endforeach()
  if(BUILD_TEST)
    set(TEST_NAME "test_${name}")
    add_executable(${TEST_NAME} "${TEST_NAME}.cc")
    add_test(${TEST_NAME} ${TEST_NAME})
    target_link_libraries(${TEST_NAME} polytope)
    set_tests_properties(${TEST_NAME} 
      PROPERTIES PASS_REGULAR_EXPRESSION PASS
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests)
  endif()
endmacro()