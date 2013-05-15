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
  for (int bid = 0; bid < 10; ++bid){
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

   // {
   //   double pts[20] = {1.461680, 0.154941,
   //                     0.200000, 1.700000,
   //                     1.100000, 0.496183,
   //                     0.300000, 0.138247,
   //                     1.955810, 1.878520,
   //                     0.150000, 0.850000,
   //                     1.300000, 0.909823,
   //                     0.571635, 1.057540,
   //                     0.143394, 1.422660,
   //                     1.850000, 0.950000};
   //   vector<double> points(20);
   //   for (unsigned i = 0; i != 20; ++i)  points[i] = pts[i];
   //   TriangleTessellator<double> tessellator;
   //   Boundary2D<double> boundary;
   //   boundary.setDefaultBoundary(3);
   //   cout << "points = [" << endl;
   //   for (int i = 0; i != points.size()/2; ++i){
   //      cout << "[" << points[2*i  ] 
   //           << "," << points[2*i+1] << "],";
   //   }
   //   cout << endl << "]" << endl << endl;
   //   Tessellation<2,double> mesh; 
   //   // The unbounded tessellation
   //   tessellator.tessellate(points, mesh);
   //   cout << "ucells = [";
   //   for (unsigned i = 0; i != mesh.cells.size(); ++i) {
   //     cout << "[" << endl;
   //     for (vector<int>::const_iterator faceItr = mesh.cells[i].begin();
   //          faceItr != mesh.cells[i].end(); ++faceItr) {
   //       const unsigned iface = *faceItr < 0 ? ~(*faceItr) : *faceItr;
   //       POLY_ASSERT(iface < mesh.faceCells.size());
   //       const unsigned inode = *faceItr < 0 ? mesh.faces[iface][0] : mesh.faces[iface][1];
   //       POLY_ASSERT(inode < mesh.nodes.size());
   //       cout << "[" << mesh.nodes[2*inode] << "," << mesh.nodes[2*inode+1] << "],";
   //     }
   //     cout << endl << "],";
   //   }
   //   cout << "]" << endl;
   //   // The bounded tessellation
   //   mesh.clear();
   //   tessellator.tessellate(points, boundary.mPLCpoints, boundary.mPLC, mesh);
   //   cout << "bcells = [";
   //   for (unsigned i = 0; i != mesh.cells.size(); ++i) {
   //     cout << "[" << endl;
   //     for (vector<int>::const_iterator faceItr = mesh.cells[i].begin();
   //          faceItr != mesh.cells[i].end(); ++faceItr) {
   //       const unsigned iface = *faceItr < 0 ? ~(*faceItr) : *faceItr;
   //       POLY_ASSERT(iface < mesh.faceCells.size());
   //       const unsigned inode = *faceItr < 0 ? mesh.faces[iface][0] : mesh.faces[iface][1];
   //       POLY_ASSERT(inode < mesh.nodes.size());
   //       cout << "[" << mesh.nodes[2*inode] << "," << mesh.nodes[2*inode+1] << "],";
   //     }
   //     cout << endl << "],";
   //   }
   //   cout << "]" << endl;
   // }


   {
     cout << "\nTriangle Tessellator:\n" << endl;
     TriangleTessellator<double> tessellator;
     testAllBoundaries(tessellator);
     //outputResult(tessellator,0,10);
   }
   
#if HAVE_BOOST_VORONOI
   {
     cout << "\nBoost Tessellator:\n" << endl;
     BoostTessellator<double> tessellator;
     testAllBoundaries(tessellator);
     //outputResult(tessellator,8,6);
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
