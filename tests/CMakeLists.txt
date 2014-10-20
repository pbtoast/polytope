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


# The POLYTOPE_ADD_TEST macro
include(${PROJECT_SOURCE_DIR}/cmake/Modules/PolytopeAddTests.cmake)


#-----------------------------------------------------------------------------
# Polytope Utility Tests
#-----------------------------------------------------------------------------
POLYTOPE_ADD_TEST( "serialize"       "" )
POLYTOPE_ADD_TEST( "convexHull_2d"   "" )
POLYTOPE_ADD_TEST( "convexHull_3d"   "" )
POLYTOPE_ADD_TEST( "nearestPoint_2d" "" )
POLYTOPE_ADD_TEST( "within_2d"       "" )
POLYTOPE_ADD_TEST( "intersect_2d"    "" )
POLYTOPE_ADD_TEST( "convexIntersect" "" )
POLYTOPE_ADD_TEST( "PLC_CSG"         "" )
POLYTOPE_ADD_TEST( "hash"            "" )

#-----------------------------------------------------------------------------
# Serial Tessellator Tests
#-----------------------------------------------------------------------------

# Voro++ Tessellator tests
#
# NOTE: Temporarily suspending their build until more work can be
#       done making their output toplogically consistent
#
#POLYTOPE_ADD_TEST( "VoroPP_2d"                   ""              )
#POLYTOPE_ADD_TEST( "VoroPP_3d"                   ""              )


# 2D Triangle and Boost Tessellator tests
POLYTOPE_ADD_TEST( "UnitSquare"                  "BOOST"         )
POLYTOPE_ADD_TEST( "MeshEditor"                  "BOOST"         )
POLYTOPE_ADD_TEST( "Degenerate"                  "TRIANGLE"      )
POLYTOPE_ADD_TEST( "BoundaryJitter"              "TRIANGLE"      )
POLYTOPE_ADD_TEST( "RandomPoints"                "TRIANGLE"      )
POLYTOPE_ADD_TEST( "Area"                        "TRIANGLE"      )
POLYTOPE_ADD_TEST( "OrphanedCell"                "TRIANGLE"      )
POLYTOPE_ADD_TEST( "TwoGenerators"               "TRIANGLE"      )
POLYTOPE_ADD_TEST( "Centroidal"                  "TRIANGLE"      )
POLYTOPE_ADD_TEST( "Unbounded"                   "TRIANGLE"      )
POLYTOPE_ADD_TEST( "UnboundedToBounded"          "TRIANGLE"      )
POLYTOPE_ADD_TEST( "OrphanCases"                 "TRIANGLE"      )
POLYTOPE_ADD_TEST( "StarBoundary"                "TRIANGLE"      )
POLYTOPE_ADD_TEST( "Circumcenters"               "TRIANGLE"      )
POLYTOPE_ADD_TEST( "BoostTessellator"            "BOOST_VORONOI" )
POLYTOPE_ADD_TEST( "RotationTests"               "TRIANGLE"      )
POLYTOPE_ADD_TEST( "SolidRotationAroundHoles"    "TRIANGLE"      )
POLYTOPE_ADD_TEST( "TiltedLattice"               "TRIANGLE"      )
#POLYTOPE_ADD_TEST( "plot"                        "BOOST_VORONOI" )

# 3D Tetgen Tessellator tests
POLYTOPE_ADD_TEST( "TetgenTessellator"           "TETGEN"        )



#-----------------------------------------------------------------------------
# Distributed Tessellator Tests
#-----------------------------------------------------------------------------
POLYTOPE_ADD_DISTRIBUTED_TEST( "DistributedUnitSquare"     "TRIANGLE" "4")
POLYTOPE_ADD_DISTRIBUTED_TEST( "SerialDistributedTriangle" "TRIANGLE" "4")
POLYTOPE_ADD_DISTRIBUTED_TEST( "FailedCommunication"       "TRIANGLE" "4")
POLYTOPE_ADD_DISTRIBUTED_TEST( "DistributedRandomPoints"   "TRIANGLE" "4")
POLYTOPE_ADD_DISTRIBUTED_TEST( "DistributedRotationTests"  "TRIANGLE" "4")
POLYTOPE_ADD_DISTRIBUTED_TEST( "DistributedCentroidal"     "TRIANGLE" "4")
POLYTOPE_ADD_DISTRIBUTED_TEST( "DistributedVoroPP_2d"      ""         "4")
POLYTOPE_ADD_DISTRIBUTED_TEST( "DistributedVoroPP_3d"      ""         "4")
