// test_Centroidal
//
// Perform repeated Lloyd iterations on the points inside a specified boundary
// to converge (slowly) to a centroidal Voronoi

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
// lloydTest
// -----------------------------------------------------------------------
void lloydTest(Tessellator<2,double>& tessellator) {
  const unsigned nPoints = 1000;     // Number of generators
  const unsigned nIter   = 100;      // Number of iterations
  const int btype = 2;
  
  string testName = "Centroidal_LloydTest_" + tessellator.name();

  // Set up boundary and disperse random generator locations
  Boundary2D<double> boundary;
  boundary.setDefaultBoundary(btype);
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
  tessellator.tessellate(points, boundary.mPLCpoints, boundary.mPLC, mesh);

  unsigned iter = 0;
  outputMesh(mesh, testName, points, iter);
  while (iter != nIter) {
    lloyd(mesh,points);
    ++iter;
    mesh.clear();
    tessellator.tessellate(points, boundary.mPLCpoints, boundary.mPLC, mesh);
    outputMesh(mesh, testName, points, iter);
  }
}

// -----------------------------------------------------------------------
// cleaningTest
// -----------------------------------------------------------------------
void cleaningTest(Tessellator<2,double>& tessellator) {
  const unsigned nPoints = 100;     // Number of generators
  const unsigned nIter   = 100;     // Number of iterations
  const double edgeTol = 0.001;     // Relative small-edge tolerance
  const int btype = 3;

  string testName = "Centroidal_CleaningTest_" + tessellator.name();

  // Set up boundary and disperse random generator locations
  Boundary2D<double> boundary;
  boundary.setDefaultBoundary(btype);  
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
  MeshEditor<2, double> meshEditor(mesh);
  tessellator.tessellate(points, boundary.mPLCpoints, boundary.mPLC, mesh);
  
  unsigned iter = 0;
  outputMesh(mesh, testName, points, iter);
  while (iter != nIter) {
    meshEditor.cleanEdges(edgeTol);
    lloyd(mesh,points);
    // meshEditor.cleanEdges(edgeTol);
    ++iter;
    mesh.clear();
    tessellator.tessellate(points, boundary.mPLCpoints, boundary.mPLC, mesh);
    outputMesh(mesh, testName, points, iter);
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
    lloydTest(tessellator);
    cleaningTest(tessellator);
  }
#endif   

#if HAVE_BOOST_VORONOI
  {
    cout << "\nBoost Tessellator:\n" << endl;
    BoostTessellator<double> tessellator;
    lloydTest(tessellator);
    cleaningTest(tessellator);
  }
#endif

  cout << "PASS" << endl;
   
#if HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
