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


# The polytope_add_test macro
include(PolytopeAddTests)

#-----------------------------------------------------------------------------
# Polytope Utility Tests
#-----------------------------------------------------------------------------
polytope_add_test( "serialize"       "" )
polytope_add_test( "convexHull_2d"   "" )
polytope_add_test( "convexHull_3d"   "" )
polytope_add_test( "nearestPoint_2d" "" )
polytope_add_test( "within_2d"       "" )
polytope_add_test( "intersect_2d"    "" )
polytope_add_test( "convexIntersect" "" )
#polytope_add_test( "PLC_CSG"         "" )  #  <-- Suspending 'til we think it works
polytope_add_test( "hash"            "" )

#-----------------------------------------------------------------------------
# Serial Tessellator Tests
#-----------------------------------------------------------------------------

# Voro++ Tessellator tests
#
# NOTE: Temporarily suspending their build until more work can be
#       done making their output toplogically consistent
#
#polytope_add_test( "VoroPP_2d"                   ""              )
#polytope_add_test( "VoroPP_3d"                   ""              )

# 2D Triangle and Boost Tessellator tests
POLYTOPE_ADD_TEST( "UnitSquare"                  ""              )
POLYTOPE_ADD_TEST( "MeshEditor"                 ""               )
POLYTOPE_ADD_TEST( "Degenerate"                  ""              )
POLYTOPE_ADD_TEST( "BoundaryJitter"              ""              )
POLYTOPE_ADD_TEST( "RandomPoints"                ""              )
POLYTOPE_ADD_TEST( "Area"                        ""              )
POLYTOPE_ADD_TEST( "OrphanedCell"                ""              )
POLYTOPE_ADD_TEST( "TwoGenerators"               ""              )
POLYTOPE_ADD_TEST( "Centroidal"                  ""              )
POLYTOPE_ADD_TEST( "Unbounded"                   ""              )
POLYTOPE_ADD_TEST( "UnboundedToBounded"          ""              )
POLYTOPE_ADD_TEST( "OrphanCases"                 ""              )
POLYTOPE_ADD_TEST( "StarBoundary"                ""              )
POLYTOPE_ADD_TEST( "Circumcenters"               ""              )
POLYTOPE_ADD_TEST( "BoostTessellator"            "BOOST_VORONOI" )
POLYTOPE_ADD_TEST( "RotationTests"               ""              )
POLYTOPE_ADD_TEST( "SolidRotationAroundHoles"    ""              )
POLYTOPE_ADD_TEST( "TiltedLattice"               ""              )
POLYTOPE_ADD_TEST( "ProjectionIntersection"      ""              )
POLYTOPE_ADD_TEST( "BoostGeometry"               "BOOST"         )
#POLYTOPE_ADD_TEST( "plot"                        "TRIANGLE"      )
#POLYTOPE_ADD_TEST( "AspectRatio"                 "TRIANGLE"      )

# 3D Tetgen Tessellator tests
polytope_add_test( "TetgenTessellator"           "TETGEN"        )

# Distributed Tessellator Tests
#-----------------------------------------------------------------------------
POLYTOPE_ADD_DISTRIBUTED_TEST( "DistributedUnitSquare"     ""         "4")
POLYTOPE_ADD_DISTRIBUTED_TEST( "SerialDistributedTriangle" "TRIANGLE" "4")
POLYTOPE_ADD_DISTRIBUTED_TEST( "FailedCommunication"       "TRIANGLE" "4")
POLYTOPE_ADD_DISTRIBUTED_TEST( "DistributedRandomPoints"   ""         "4")
POLYTOPE_ADD_DISTRIBUTED_TEST( "DistributedRotationTests"  ""         "4")
POLYTOPE_ADD_DISTRIBUTED_TEST( "DistributedCentroidal"     ""         "4")
POLYTOPE_ADD_DISTRIBUTED_TEST( "DistributedCollinear"      "TRIANGLE" "4")
#POLYTOPE_ADD_DISTRIBUTED_TEST( "DistributedVoroPP_2d"      ""         "4")
#POLYTOPE_ADD_DISTRIBUTED_TEST( "DistributedVoroPP_3d"      ""         "4")
