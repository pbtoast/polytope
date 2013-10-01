#-----------------------------------------------------------------------------
# Adding a C test to polytope_c
# 
# This test directory assumes all tests follow the naming convention
#      "test_<name>.c"
#
# To add a test, the user must register the test <name> and a list of the
# dependencies needed to build it: 
#      POLYTOPE_C_ADD_TEST( name dependency_list )
#
# The current set of dependencies is
#      MPI           : MPI is used
#      TETGEN        : Tetgen Tessellator is used
#      BOOST         : Boost library is used
#      TRIANGLE      : Triangle Tessellator is used
#      BOOST_VORONOI : Boost v1.52 or greater is used having the
#                      Boost.Polygon Voronoi library
#
# NOTES:
# (1) An empty string in the dependency list means the test will always build
# (2) If you need Boost and Triangle, you need only put "TRIANGLE"
# (3) If a test has multiple dependencies, submit the dependency list
#     as a single string, with individual entries separated by a semi-colon
#     (EXAMPLE: if both TETGEN and MPI are needed, the dependency list
#               should read "TETGEN;MPI"
#-----------------------------------------------------------------------------


# The POLYTOPE_C_ADD_TEST macro
include(${PROJECT_SOURCE_DIR}/cmake/Modules/PolytopeCAddTests.cmake)

# 2D Triangle and Boost Tessellator tests
POLYTOPE_C_ADD_TEST( "boost_tessellator"            "BOOST_VORONOI" )

# 3D Tetgen Tessellator tests
#POLYTOPE_C_ADD_TEST( "tetgen_tessellator"           "TETGEN"        )

#-----------------------------------------------------------------------------
# Distributed Tessellator Tests
#-----------------------------------------------------------------------------
#POLYTOPE_C_ADD_DISTRIBUTED_TEST( "distributed_unit_square"     "TRIANGLE" "4")
#POLYTOPE_C_ADD_DISTRIBUTED_TEST( "serial_distributed_triangle" "TRIANGLE" "4")
#POLYTOPE_C_ADD_DISTRIBUTED_TEST( "failed_communication"        "TRIANGLE" "4")
#POLYTOPE_C_ADD_DISTRIBUTED_TEST( "distributed_random_points"   "TRIANGLE" "4")
#POLYTOPE_C_ADD_DISTRIBUTED_TEST( "distributed_rotation_tests"  "TRIANGLE" "4")
#POLYTOPE_C_ADD_DISTRIBUTED_TEST( "distributed_centroidal"      "TRIANGLE" "4")
