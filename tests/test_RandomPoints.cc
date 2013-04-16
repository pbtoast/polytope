// test_RandomPoints
//
// Stress test for meshing complicated PLC boundaries with/without holes.
// Iterate over each of the default boundaries defined in Boundary2D.hh
// and tessellate using N randomly-distributed generators for N=10,100,1000.
// Can test both Triangle and Voro++ 2D tessellators. Voro++ has been
// commented out since it currently lacks PLC capabilities.

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
// outputResult
// -----------------------------------------------------------------------
void outputResult(Tessellator<2,double>& tessellator,
		  const int bType,
		  const unsigned nPoints) {
  // output name
  ostringstream os;
  os << "RandomPoints_" << tessellator.name() << "_" << bType;
  string testName = os.str();

  // Boundary data
  Boundary2D<double> boundary;
  boundary.setDefaultBoundary(bType);
  
  // Generator data
  Generators<2,double> generators( boundary );
  generators.randomPoints(nPoints);
  
  Tessellation<2,double> mesh;
  tessellate2D(generators.mPoints,boundary,tessellator,mesh);
  POLY_ASSERT( mesh.cells.size() == nPoints );
  double area = computeTessellationArea(mesh);
  cout << "Tessellation Area = " << area << endl;
  cout << "Relative error    = " << (boundary.mArea-area)/boundary.mArea << endl;  
  outputMesh(mesh, testName, generators.mPoints);
}


// -----------------------------------------------------------------------
// testBoundary
// -----------------------------------------------------------------------
void testBoundary(Boundary2D<double>& boundary,
                  Tessellator<2,double>& tessellator) {
  // output name
  ostringstream os;
  os << "RandomPoints_" << tessellator.name() << "_" << boundary.mType;
  string testName = os.str();

  Generators<2,double> generators( boundary );
  unsigned nPoints = 1;
  Tessellation<2,double> mesh;
  for( unsigned n = 0; n < 3; ++n ){
    POLY_ASSERT(mesh.empty());
    nPoints = nPoints * 10;
    
    cout << nPoints << " points...";
    generators.randomPoints( nPoints );      
    tessellate2D(generators.mPoints,boundary,tessellator,mesh);
    cout << "got meshed!" << endl;
    
    outputMesh(mesh, testName, generators.mPoints, n);
    mesh.clear();   
  }
}

// -----------------------------------------------------------------------
// testAllBoundaries
// -----------------------------------------------------------------------
void testAllBoundaries(Tessellator<2,double>& tessellator) {
  for (int bid = 0; bid < 9; ++bid){
    cout << "Testing boundary type " << bid << endl;
    Boundary2D<double> boundary;
    boundary.setDefaultBoundary(bid);
    testBoundary(boundary, tessellator);
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

   {
     cout << "\nTriangle Tessellator:\n" << endl;
     TriangleTessellator<double> tessellator;
     testAllBoundaries(tessellator);
     //outputResult(tessellator,3,20);
   }
   
#if HAVE_BOOST_VORONOI
   {
     cout << "\nBoost Tessellator:\n" << endl;
     BoostTessellator<double> tessellator;
     testAllBoundaries(tessellator);
   }
#endif      

   // NOTE: Voro++ currently lacks PLC boundary capabilities
   //
   // cout << "\nVoro 2D Tessellator:\n" << endl;
   // VoroPP_2d<double> voro;
   // testAllBoundaries(voro);

   cout << "PASS" << endl;

#if HAVE_MPI
   MPI_Finalize();
#endif
   return 0;
}
