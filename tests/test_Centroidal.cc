
#include <iostream>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>
#include <sstream>
#include <algorithm>

#include "polytope.hh"
#include "Boundary2D.hh"
#include "Generators.hh"
#include "polytope_test_utilities.hh"
#include "polytope_geometric_utilities.hh"

#if HAVE_MPI
#include "mpi.h"
#endif

using namespace std;
using namespace polytope;

// -----------------------------------------------------------------------
// lloyd
// -----------------------------------------------------------------------
void lloyd(Tessellation<2,double>& mesh,
           vector<double>& points) {
   for (unsigned i = 0; i < mesh.cells.size(); ++i){
      double cent[2], area;
      geometry::computeCellCentroidAndSignedArea(mesh, i, 1.0e-12, cent, area);
      points[2*i  ] = 0.5*(points[2*i  ] + cent[0]);
      points[2*i+1] = 0.5*(points[2*i+1] + cent[1]);
   }
}


// -----------------------------------------------------------------------
// test
// -----------------------------------------------------------------------
void test(Tessellator<2,double>& tessellator) {
  const unsigned nPoints = 100;     // Number of generators
  const unsigned nIter   = 100;     // Number of iterations

  string testName = "Centroidal_" + tessellator.name();

  // Set up boundary and disperse random generator locations
  Boundary2D<double> boundary;
  boundary.setDonut();
  Generators<2,double> generators(boundary);
  generators.randomPoints(nPoints);
  std::vector<double> points;
  for (unsigned i = 0; i != nPoints; ++i) {
    if (boundary.testInside(&generators.mPoints[2*i])) {
      std::copy(&generators.mPoints[2*i], &generators.mPoints[2*i+2], std::back_inserter(points));
    }
  }
   
  // Initialize mesh and tessellator
  Tessellation<2,double> mesh;

  for( int iter = 0; iter < nIter; ++iter ){
    mesh.clear();
    tessellator.tessellate(points, boundary.mPLCpoints, boundary.mPLC, mesh);
    outputMesh(mesh, testName, points, iter);
    lloyd(mesh,points);
  }
}


// -----------------------------------------------------------------------
// main
// -----------------------------------------------------------------------
int main(int argc, char** argv)
{
#if HAVE_MPI
  MPI_Init(&argc, &argv);
#endif
   
#if HAVE_TRIANGLE
  {
    cout << "\nTriangle Tessellator:\n" << endl;
    TriangleTessellator<double> tessellator;
    test(tessellator);
  }
#endif   

#if HAVE_BOOST_VORONOI
  {
    cout << "\nBoost Tessellator:\n" << endl;
    BoostTessellator<double> tessellator;
    test(tessellator);
  }
#endif

  cout << "PASS" << endl;
   
#if HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
