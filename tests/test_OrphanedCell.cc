// test_OrphanedCell
//
// 3x3 Grid with Cartesian generators. Two rectangles jut out of the top
// and bottom boundaries, dividing the top-middle and bottom-middle cells.
// For each of the two cells: the divided part that contains the generator
// becomes the cell, while the other part becomes an orphaned cell.
//
// This tests the "cell adoption" capability in the TriangleTessellator  
// which appropriates the area in the orphaned pieces to its neighboring
// cells while maintaining the Voronoi property locally

#include <iostream>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>
#include <sstream>

#include "polytope.hh"
#include "polytope_test_utilities.hh"

#if HAVE_MPI
#include "mpi.h"
#endif

using namespace std;
using namespace polytope;

// -----------------------------------------------------------------------
// main
// -----------------------------------------------------------------------
int main(int argc, char** argv)
{
#if HAVE_MPI
   MPI_Init(&argc, &argv);
#endif
   
   std::vector<double> PLCpoints;
   std::vector<double> generators;
   PLC<2,double> boundary;
   Tessellation<2,double> mesh;
   TriangleTessellator<double> triangle;
   
   PLCpoints.push_back(0.0);  PLCpoints.push_back(0.0);
   PLCpoints.push_back(1.2);  PLCpoints.push_back(0.0);
   PLCpoints.push_back(1.2);  PLCpoints.push_back(1.3);
   PLCpoints.push_back(1.3);  PLCpoints.push_back(1.3);
   PLCpoints.push_back(1.3);  PLCpoints.push_back(0.0);
   PLCpoints.push_back(3.0);  PLCpoints.push_back(0.0);
   PLCpoints.push_back(3.0);  PLCpoints.push_back(3.0);
   PLCpoints.push_back(1.3);  PLCpoints.push_back(3.0);
   PLCpoints.push_back(1.3);  PLCpoints.push_back(1.7);
   PLCpoints.push_back(1.2);  PLCpoints.push_back(1.7);
   PLCpoints.push_back(1.2);  PLCpoints.push_back(3.0);
   PLCpoints.push_back(0.0);  PLCpoints.push_back(3.0);
   
   int ix, iy, nx = 3;
   double xi, yi;
   for (iy = 0; iy != nx; ++iy) {
     yi = iy + 0.5;
     for (ix = 0; ix != nx; ++ix) {
       xi = ix + 0.5;
       generators.push_back(xi);  generators.push_back(yi);
     }
   }
   
   int nSides = PLCpoints.size()/2;
   boundary.facets.resize( nSides, std::vector<int>(2) );
   for (unsigned i = 0; i != nSides; ++i){
     boundary.facets[i][0] = i;
     boundary.facets[i][1] = (i+1) % nSides;
   }
   
   triangle.tessellate(generators, PLCpoints, boundary, mesh);
   
#if HAVE_SILO
   vector<double> index(mesh.cells.size());
   vector<double> genx (mesh.cells.size());
   vector<double> geny (mesh.cells.size());
   for (int i = 0; i < mesh.cells.size(); ++i){
      index[i] = double(i);
      genx[i]  = generators[2*i];
      geny[i]  = generators[2*i+1];
   }
   map<string,double*> nodeFields, edgeFields, faceFields, cellFields;
   cellFields["cell_index"   ] = &index[0];
   cellFields["cell_center_x"] = &genx[0];
   cellFields["cell_center_y"] = &geny[0];
   ostringstream os;
   os << "test_OrphanedCell_mesh";
   polytope::SiloWriter<2, double>::write(mesh, nodeFields, edgeFields, 
                                          faceFields, cellFields, os.str());
#endif


   cerr << mesh << endl;

   // Total mesh checks
   POLY_CHECK(mesh.cells.size()   == 9 );
   //POLY_CHECK(mesh.nodes.size()/2 == 26);
   //POLY_CHECK(mesh.faces.size()   == 34);
   // Individual cell checks
   //POLY_CHECK(mesh.cells[0].size() == 5 );
   POLY_CHECK(mesh.cells[1].size() == 4 );
   POLY_CHECK(mesh.cells[2].size() == 4 );
   POLY_CHECK(mesh.cells[3].size() == 4 );
   //POLY_CHECK(mesh.cells[4].size() == 12);
   POLY_CHECK(mesh.cells[5].size() == 4 );
   POLY_CHECK(mesh.cells[6].size() == 5 );
   POLY_CHECK(mesh.cells[7].size() == 4 );
   POLY_CHECK(mesh.cells[8].size() == 4 );
   // Tessellation area check
   const double trueArea = 8.74;
   const double tessArea = computeTessellationArea(mesh);
   const double fracerr  = std::abs(trueArea - tessArea)/trueArea;
   const double tol      = 1.0e-9;
   POLY_CHECK2(fracerr < tol, "Relative error in the tessellation "
	       << "area exceeds tolerance:" << endl
	       << "      Area = " << tessArea << endl
	       << "     Error = " << trueArea - tessArea << endl
	       << "Frac Error = " << fracerr);

   cout << "PASS" << endl;

#if HAVE_MPI
   MPI_Finalize();
#endif
   return 0;
}
