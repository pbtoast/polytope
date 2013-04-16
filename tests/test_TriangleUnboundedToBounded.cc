// test_TriangleUnboundedToBounded

#include <iostream>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>
#include <sstream>

#include "polytope.hh"
#include "Boundary2D.hh"
#include "Generators.hh"
#include "polytope_test_utilities.hh"

#if HAVE_MPI
#include "mpi.h"
#endif

using namespace std;
using namespace polytope;


// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
int main(int argc, char** argv)
{
#if HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  // unsigned nPoints = 5;
  // double pts[10] = {0.0, 0.0,
  //                   1.0, 0.0,
  //                   1.0, 1.0,
  //                   0.0, 1.0,
  //                   0.2, 0.5};

  // std::vector<double> points;
  // for (unsigned i=0; i<nPoints; ++i){
  //    points.push_back(pts[2*i  ]);
  //    points.push_back(pts[2*i+1]);
  // }

  const string testName = "TriangleUnboundedToBounded";
  const int bType = 3;
  
  // Boundary data
  Boundary2D<double> boundary;
  boundary.setDefaultBoundary(bType);

  // Generator data
  Generators<2,double> generators(boundary);
  generators.randomPoints(20);

  Tessellation<2,double> mesh;
  TriangleTessellator<double> tessellator;

  // Do the unbounded tessellation
  tessellator.tessellate( generators.mPoints, mesh );
  outputMesh(mesh, testName, generators.mPoints, 1);

  // const unsigned nSides = 4;
  // std::vector<double> PLCpoints;
  // PLCpoints.push_back(0.0);  PLCpoints.push_back(0.0);
  // PLCpoints.push_back(1.0);  PLCpoints.push_back(0.0);
  // PLCpoints.push_back(1.0);  PLCpoints.push_back(1.0);
  // PLCpoints.push_back(0.0);  PLCpoints.push_back(1.0);
  
  // PLC<2,double> geometry;
  // geometry.facets.resize( nSides, std::vector<int>(2) );
  // for (unsigned i = 0; i != nSides; ++i){
  //   geometry.facets[i][0] = i;
  //   geometry.facets[i][1] = (i+1) % nSides;
  // }
  
  // Now do the bounded version
  mesh.clear();
  tessellator.tessellate( generators.mPoints, boundary.mPLCpoints, boundary.mPLC, mesh );
  outputMesh(mesh, testName, generators.mPoints, 2);

  
  cout << "PASS" << endl;

#if HAVE_MPI
   MPI_Finalize();
#endif
   return 0;
}
