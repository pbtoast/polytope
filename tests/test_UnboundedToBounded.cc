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
// runTest
// -----------------------------------------------------------------------
void runTest(Tessellator<2,double>& tessellator) {
  ostringstream os;
  os << "UnboundedToBounded_" << tessellator.name();
  const string testName = os.str();
  
  const int bType = 3;  // M-shape with holes
  
  // Boundary data
  Boundary2D<double> boundary;
  boundary.setDefaultBoundary(bType);

  // Generator data
  Generators<2,double> generators(boundary);
  generators.randomPoints(20);

  Tessellation<2,double> mesh;
  
  // Do the unbounded tessellation
  tessellator.tessellate( generators.mPoints, mesh );
  outputMesh(mesh, testName, generators.mPoints, 1);
  
  // Now do the bounded version
  mesh.clear();
  tessellator.tessellate( generators.mPoints, boundary.mPLCpoints, boundary.mPLC, mesh );
  outputMesh(mesh, testName, generators.mPoints, 2);
}

// -----------------------------------------------------------------------
// main
// -----------------------------------------------------------------------
int main(int argc, char** argv)
{
#if HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  {
    cout << "\nTriangle Tessellator:\n" << endl;
    TriangleTessellator<double> tessellator;
    runTest(tessellator);
  }
   
#if HAVE_BOOST_VORONOI
  {
    cout << "\nBoost Tessellator:\n" << endl;
    BoostTessellator<double> tessellator;
    runTest(tessellator);
  }
#endif      
  
  cout << "PASS" << endl;

#if HAVE_MPI
  MPI_Finalize();
#endif
   return 0;
}
