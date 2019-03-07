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
polytope_add_test( "UnitSquare"                  "BOOST"         )
polytope_add_test( "MeshEditor"                  "BOOST"         )
polytope_add_test( "Degenerate"                  "TRIANGLE"      )
polytope_add_test( "BoundaryJitter"              "TRIANGLE"      )
polytope_add_test( "RandomPoints"                "TRIANGLE"      )
polytope_add_test( "Area"                        "TRIANGLE"      )
polytope_add_test( "OrphanedCell"                "TRIANGLE"      )
polytope_add_test( "TwoGenerators"               "TRIANGLE"      )
polytope_add_test( "Centroidal"                  "TRIANGLE"      )
polytope_add_test( "Unbounded"                   "TRIANGLE"      )
polytope_add_test( "UnboundedToBounded"          "TRIANGLE"      )
polytope_add_test( "OrphanCases"                 "TRIANGLE"      )
polytope_add_test( "StarBoundary"                "TRIANGLE"      )
polytope_add_test( "Circumcenters"               "TRIANGLE"      )
polytope_add_test( "BoostTessellator"            "BOOST_VORONOI" )
polytope_add_test( "RotationTests"               "TRIANGLE"      )
polytope_add_test( "SolidRotationAroundHoles"    "TRIANGLE"      )
polytope_add_test( "TiltedLattice"               "TRIANGLE"      )
polytope_add_test( "ProjectionIntersection"      "TRIANGLE"      )
polytope_add_test( "BoostGeometry"               "BOOST"         )
#polytope_add_test( "plot"                        "TRIANGLE"      )
#polytope_add_test( "AspectRatio"                 "TRIANGLE"      )

# 3D Tetgen Tessellator tests
polytope_add_test( "TetgenTessellator"           "TETGEN"        )

# Distributed Tessellator Tests
polytope_add_distributed_test( "DistributedUnitSquare"     "TRIANGLE" "4")
polytope_add_distributed_test( "SerialDistributedTriangle" "TRIANGLE" "4")
polytope_add_distributed_test( "FailedCommunication"       "TRIANGLE" "4")
polytope_add_distributed_test( "DistributedRandomPoints"   "TRIANGLE" "4")
polytope_add_distributed_test( "DistributedRotationTests"  "TRIANGLE" "4")
polytope_add_distributed_test( "DistributedCentroidal"     "TRIANGLE" "4")
polytope_add_distributed_test( "DistributedCollinear"      "TRIANGLE" "4")
#polytope_add_distributed_test( "DistributedVoroPP_2d"      ""         "4")
#polytope_add_distributed_test( "DistributedVoroPP_3d"      ""         "4")
