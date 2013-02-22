#-----------------------------------------------------------------------------
# Adding a C++ test to Polytope
# 
# The polytope test directory assumes all tests follow the naming convention
#      "test_<name>.cc"
#
# To add a test, the user must register the test <name> and a list of the
# dependencies needed to build it: 
#      POLYTOPE_ADD_TEST( name dependency_list )
#
# The current set of dependencies is
#      MPI       : MPI is used
#      TETGEN    : Tetgen Tessellator is used
#      BOOST     : Boost library is used
#      TRIANGLE  : Triangle Tessellator is used
#
# NOTES:
# (1) An empty string in the dependency list means the test will always build
# (2) If you need Boost and Triangle, you need only put "TRIANGLE"
# (3) If a test has multiple dependencies, submit the dependency list
#     as a single string, with individual entries separated by a semi-colon
#     (EXAMPLE: if both TETGEN and MPI are needed, the dependency list
#               should read "TETGEN;MPI"
#-----------------------------------------------------------------------------

# The POLYTOPE_ADD_TEST macro
include(${PROJECT_SOURCE_DIR}/cmake/Modules/PolytopeAddTests.cmake)

# Tests not involving the tessellators
POLYTOPE_ADD_TEST( "serialize"       "" )
POLYTOPE_ADD_TEST( "convexHull_2d"   "" )
POLYTOPE_ADD_TEST( "convexHull_3d"   "" )
POLYTOPE_ADD_TEST( "nearestPoint_2d" "" )
POLYTOPE_ADD_TEST( "within_2d"       "" )
POLYTOPE_ADD_TEST( "intersect_2d"    "" )

# Serial tessellator tests
POLYTOPE_ADD_TEST( "UnitSquare"          "BOOST"    )
POLYTOPE_ADD_TEST( "Degenerate"          "BOOST"    )
POLYTOPE_ADD_TEST( "TriangleTessellator" "TRIANGLE" )
POLYTOPE_ADD_TEST( "RandomPoints"        "TRIANGLE" )
POLYTOPE_ADD_TEST( "Area"                "TRIANGLE" )
POLYTOPE_ADD_TEST( "OrphanedCell"        "TRIANGLE" )
POLYTOPE_ADD_TEST( "TwoGenerators"       "TRIANGLE" )
POLYTOPE_ADD_TEST( "TetgenTessellator"   "TETGEN"   )
POLYTOPE_ADD_TEST( "VoroPP_2d"           ""         )
POLYTOPE_ADD_TEST( "VoroPP_3d"           ""         )
POLYTOPE_ADD_TEST( "Centroidal"          "TRIANGLE" )
POLYTOPE_ADD_TEST( "plcOrder"            "TRIANGLE" )
POLYTOPE_ADD_TEST( "TriangleUnbounded"   "TRIANGLE" )
POLYTOPE_ADD_TEST( "OrphanCases"         "TRIANGLE" )

# Distributed tessellator tests
POLYTOPE_ADD_TEST( "DistributedTriangle"       "MPI;TRIANGLE" )
POLYTOPE_ADD_TEST( "SerialDistributedTriangle" "MPI;TRIANGLE" )
POLYTOPE_ADD_TEST( "FailedCommunication"       "MPI;TRIANGLE" )
POLYTOPE_ADD_TEST( "DistributedRandomPoints"   "MPI;TRIANGLE" )
POLYTOPE_ADD_TEST( "DistributedVoroPP_2d"      "MPI"          )
POLYTOPE_ADD_TEST( "DistributedVoroPP_3d"      "MPI"          )
