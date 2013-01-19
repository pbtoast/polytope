// test_Area
//
// Compare the total area of the polytope tessellation to the area of the
// bounding region. Iterates over all default boundaries in Boundary2D.hh
// and computes tessellation areas for N randomly-distributed generators,
// where N = 10, 100, 1000, and 10000. Areas are computed using 
// Boost.Geometry
//
// NOTE: Area discrepancies are due to cells divided by PLC boundaries.
//       The part of the divided cell that contains the generator is
//       kept; the other part is lost and is not represented as a mesh
//       cell. Adding generators does not always rectify this, for instance
//       if the boundary is a cardioid.


#include <iostream>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>
#include <sstream>

#include "polytope.hh"
#include "polytope_test_utilities.hh"
#include "Generators.hh"

#if HAVE_MPI
// extern "C" {
#include "mpi.h"
// }
#endif

using namespace std;
using namespace polytope;

// -----------------------------------------------------------------------
// testBoundary
// -----------------------------------------------------------------------
double testBoundary(Boundary2D<double>& boundary,
                    Tessellator<2,double>& tessellator)
{
   Generators<2,double> generators( boundary );
   unsigned nPoints = 1;
   Tessellation<2,double> mesh;
   cout << "Area of boundary = " << boundary.mArea << endl;
   double result = 0.0;
   for( unsigned n = 0; n < 4; ++n ){
      POLY_ASSERT( mesh.empty() );
      nPoints = nPoints * 10;
      cout << nPoints << " points..." << endl;

      generators.randomPoints( nPoints );
      tessellate2D(generators.mPoints, boundary, tessellator, mesh);

      const double area = computeTessellationArea( mesh );
      const double fracerr = std::abs(boundary.mArea - area)/boundary.mArea;
      result = std::max(result, fracerr);
      cout << "              Area  = " << area << endl;
      cout << "              Error = " << boundary.mArea - area << endl;
      cout << "   Fractional error = " << fracerr << endl;

      mesh.clear();   
   }
   return result;
}


// -----------------------------------------------------------------------
// testAllBoundaries
// -----------------------------------------------------------------------
double testAllBoundaries(Tessellator<2,double>& tessellator)
{
   double result = 0.0;
   for (int bid = 0; bid < 7; ++bid){
      cout << endl << "Testing boundary type " << bid << endl;
      Boundary2D<double> boundary;
      boundary.computeDefaultBoundary(bid);
      result = std::max(result, testBoundary( boundary, tessellator ));
   }
   return result;
}


// -----------------------------------------------------------------------
// the test
// -----------------------------------------------------------------------
int main(int argc, char** argv)
{
#if HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  cout << "\nTriangle Tessellator:\n" << endl;
  TriangleTessellator<double> triangle;
  const double maxError = testAllBoundaries(triangle);
   
  const double tol = 0.1;
  if (maxError > tol) {
     cout << "FAIL" << endl;
  } else {
     cout << "PASS" << endl;
  }

#if HAVE_MPI
  MPI_Finalize();
#endif
   return 0;
}




// NOTE: Old implementation for computing tessellation areas. Current
//       implementation based on Boost.Geometry found in polytope_test_utilities.hh
//
// // -------------------------------------------------------------------- //
// double cellArea(std::vector<double> x, std::vector<double> y )
// {
//    double  area=0.0;
//    POLY_ASSERT( x.size() == y.size() );
//    int j=x.size()-1;
//    for (int i = 0; i < x.size(); ++i) {
//       area -= (x[j] + x[i]) * (y[j] - y[i]); 
//       j=i; }
//    return 0.5*area; 
// }
// // -------------------------------------------------------------------- //
// double tessellationArea( Tessellation<2,double>& mesh ){
//    std::vector<double> x,y;
//    double area = 0.0;
//    for (unsigned i = 0; i < mesh.cells.size(); ++i )
//    {
//       // cout << endl << "Cell " << i << " has node positions" << endl;
//       for (std::vector<int>::const_iterator faceItr = 
//               mesh.cells[i].begin(); faceItr != 
//               mesh.cells[i].end();  ++faceItr)
//       {
//          const unsigned iface = *faceItr < 0 ? ~(*faceItr) : *faceItr;
//          POLY_ASSERT( mesh.faces[iface].size() == 2 );
//          const unsigned inode = *faceItr < 0 ? mesh.faces[iface][1] :
//             mesh.faces[iface][0];
//          x.push_back( mesh.nodes[2*inode  ] );
//          y.push_back( mesh.nodes[2*inode+1] );
//       }
//       area += cellArea( x, y );
//       x.clear(); y.clear();
//    }
//    return area;
// }
// // -------------------------------------------------------------------- //
