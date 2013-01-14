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
// extern "C" {
#include "mpi.h"
// }
#endif

using namespace std;
using namespace polytope;

// -----------------------------------------------------------------------
// outputResult
// -----------------------------------------------------------------------
void outputResult(Tessellator<2,double>& tessellator,
		  int bType,
		  unsigned nPoints)
{
   Boundary2D<double> boundary;

   boundary.computeDefaultBoundary(bType);
   Generators<2,double> generators( boundary );

   generators.randomPoints( nPoints );

   Tessellation<2,double> mesh;
   tessellate2D(generators.mPoints,boundary,tessellator,mesh);
   POLY_ASSERT( mesh.cells.size() == nPoints );
   
#if HAVE_SILO
   vector<double> index( mesh.cells.size());
   for (int i = 0; i < mesh.cells.size(); ++i) index[i] = double(i);
   map<string,double*> nodeFields, edgeFields, faceFields, cellFields;
   cellFields["cell_index"] = &index[0];
   ostringstream os;
   os << "test_RandomPoints_boundary_" << bType << "_" << nPoints << "points";
   polytope::SiloWriter<2, double>::write(mesh, nodeFields, edgeFields, 
                                          faceFields, cellFields, os.str());
#endif
}


// -----------------------------------------------------------------------
// testBoundary
// -----------------------------------------------------------------------
void testBoundary(Boundary2D<double>& boundary,
                  Tessellator<2,double>& tessellator)
{
   Generators<2,double> generators( boundary );
   unsigned nPoints = 1;
   Tessellation<2,double> mesh;
   for( unsigned n = 0; n < 3; ++n ){
      POLY_ASSERT( mesh.empty() );
      nPoints = nPoints * 10;
      cout << nPoints << " points...";

      generators.randomPoints( nPoints );      
      tessellate2D(generators.mPoints,boundary,tessellator,mesh);

      POLY_ASSERT( mesh.cells.size() == nPoints );
      cout << "PASS" << endl;
#if HAVE_SILO
      vector<double> index( mesh.cells.size());
      for (int i = 0; i < mesh.cells.size(); ++i) index[i] = double(i);
      map<string,double*> nodeFields, edgeFields, faceFields, cellFields;
      cellFields["cell_index"] = &index[0];
      ostringstream os;
      os << "test_RandomPoints_boundary_" << boundary.mType << "_" << nPoints << "_points";
      polytope::SiloWriter<2, double>::write(mesh, nodeFields, edgeFields, 
                                             faceFields, cellFields, os.str());
#endif
      mesh.clear();   
   }
}

// -----------------------------------------------------------------------
// testAllBoundaries
// -----------------------------------------------------------------------
void testAllBoundaries(Tessellator<2,double>& tessellator)
{
   for (int bid = 0; bid < 7; ++bid){
      cout << "Testing boundary type " << bid << endl;
      Boundary2D<double> boundary;
      boundary.computeDefaultBoundary(bid);
      testBoundary( boundary, tessellator );
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

   cout << "\nTriangle Tessellator:\n" << endl;
   TriangleTessellator<double> triangle;
   testAllBoundaries(triangle);
   
   // NOTE: Voro++ currently lacks PLC boundary capabilities
   //
   // cout << "\nVoro 2D Tessellator:\n" << endl;
   // VoroPP_2d<double> voro;
   // testAllBoundaries(voro);

#if HAVE_MPI
   MPI_Finalize();
#endif
   return 0;
}
