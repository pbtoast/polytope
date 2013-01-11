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

using namespace std;
using namespace polytope;

// -------------------------------------------------------------------- //
template <class T>
inline std::string to_string(const T& t){
   std::stringstream ss;
   ss << t;
   return ss.str();
}
// //
// double computeArea( Tessellation<2,double>& mesh )
// {
//    double area = 0.0;
//    std::vector<std::set<unsigned> > cellToNodes = mesh.computeCellToNodes();
//    for (unsigned icell = mesh.cells.begin(); 
//         icell != mesh.cells.end();
//         ++icell){
//       std::vector<double> nodeCell;
//       for (std::set<unsigned>::const_iterator nodeItr = cellToNodes[icell].begin();
//            nodeItr != cellToNodes[icell].end(); ++nodeItr){
//          unsigned inode = *nodeItr;
//          nodeCell.push_back( mesh.nodes[2*inode  ] );
//          nodeCell.push_back( mesh.nodes[2*inode+1] );
//       }
//       area += cellArea( nodeCell );
//    }
//    return area;
// }
// // -------------------------------------------------------------------- //
double cellArea(std::vector<double> x, std::vector<double> y )
{
   double  area=0.0;
   POLY_ASSERT( x.size() == y.size() );
   int j=x.size()-1;
   for (int i = 0; i < x.size(); ++i) {
      area -= (x[j] + x[i]) * (y[j] - y[i]); 
      j=i; }
   return 0.5*area; 
}
// -------------------------------------------------------------------- //
double tessellationArea( Tessellation<2,double>& mesh ){
   std::vector<double> x,y;
   double area = 0.0;
   for (unsigned i = 0; i < mesh.cells.size(); ++i )
   {
      // cout << endl << "Cell " << i << " has node positions" << endl;
      for (std::vector<int>::const_iterator faceItr = 
              mesh.cells[i].begin(); faceItr != 
              mesh.cells[i].end();  ++faceItr)
      {
         const unsigned iface = *faceItr < 0 ? ~(*faceItr) : *faceItr;
         POLY_ASSERT( mesh.faces[iface].size() == 2 );
      
         const unsigned inode = *faceItr < 0 ? mesh.faces[iface][1] :
            mesh.faces[iface][0];
         x.push_back( mesh.nodes[2*inode  ] );
         y.push_back( mesh.nodes[2*inode+1] );
      }
      area += cellArea( x, y );
      x.clear(); y.clear();
   }
   return area;
}
// -------------------------------------------------------------------- //



// -----------------------------------------------------------------------
// testBoundary
// -----------------------------------------------------------------------
void testBoundary(Boundary2D<double>& boundary,
                  Tessellator<2,double>& tessellator)
{
   Generators<2,double> generators( boundary );
   unsigned nPoints = 1;
   Tessellation<2,double> mesh;
   cout << "Area of boundary = " << boundary.mArea << endl;
   for( unsigned n = 0; n < 4; ++n ){
      POLY_ASSERT( mesh.empty() );
      nPoints = nPoints * 10;
      cout << nPoints << " points...";

      generators.randomPoints( nPoints );
      
      tessellate2D(generators.mPoints, boundary, tessellator, mesh);

      double area = tessellationArea( mesh );
      cout << "Area of tessellation = " << area << endl;

// #if HAVE_SILO
//       std::string name = "test_RandomPoints_boundary" 
//          + to_string(boundary.mType) + "_" + to_string(nPoints) + "points";
//       vector<double> index( mesh.cells.size());
//       for (int i = 0; i < mesh.cells.size(); ++i) index[i] = double(i);
//       map<string,double*> nodeFields, edgeFields, faceFields, cellFields;
//       cellFields["cell_index"] = &index[0];
//       polytope::SiloWriter<2, double>::write(mesh, nodeFields, edgeFields, 
//                                              faceFields, cellFields, name);
// #endif
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

   
#if HAVE_MPI
   MPI_Finalize();
#endif
   return 0;
}
