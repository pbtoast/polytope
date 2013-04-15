// test_UnitSquare
//
// Mesh a unit square with (nx-by-nx) Cartesian generators for nx in [2,100].
// Perform checks on the resulting tessellation to see if it is indeed Cartesian.
// Both Triangle and Voro++ 2D tessellators are tested here.

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
// checkCartesianMesh
// -----------------------------------------------------------------------
void checkCartesianMesh(Tessellation<2,double>& mesh, unsigned nx, unsigned ny)
{
   POLY_CHECK(mesh.nodes.size()/2 == (nx + 1)*(ny + 1) );
   POLY_CHECK(mesh.cells.size()   == nx*ny );
   POLY_CHECK(mesh.faces.size()   == nx*(ny + 1) + ny*(nx + 1) );
   for (unsigned i = 0; i != nx*ny; ++i) POLY_CHECK(mesh.cells[i].size() == 4);
   
   std::vector<std::set<unsigned> > nodeCells = mesh.computeNodeCells();
   for (unsigned i = 0; i != (nx+1)*(ny+1); ++i)
   {
      POLY_CHECK( (nodeCells[i].size() == 4) ||
                  (nodeCells[i].size() == 2) ||
                  (nodeCells[i].size() == 1) );
   }
}

// -----------------------------------------------------------------------
// generateMesh
// -----------------------------------------------------------------------
void generateMesh(Tessellator<2,double>& tessellator)
{
   // Set the boundary
   Boundary2D<double> boundary;
   boundary.setUnitSquare();
   Generators<2,double> generators( boundary );
   
   for (unsigned nx = 2; nx != 100; ++nx){
      cout << "Testing nx=" << nx << endl;

      // Create generators
      std::vector<unsigned> nxny(2,nx);
      generators.cartesianPoints( nxny );
      Tessellation<2,double> mesh;
      tessellate2D(generators.mPoints,boundary,tessellator,mesh);
      
      // CHECKS:
      cout << "   num mesh nodes : " << mesh.nodes.size()/2 << endl;
      cout << "   num mesh cells : " << mesh.cells.size()   << endl;
      cout << "   num mesh faces : " << mesh.faces.size()   << endl;
      checkCartesianMesh(mesh,nx,nx);
   }
}


// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
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
    generateMesh(tessellator);
  }
#endif   

#if HAVE_BOOST_VORONOI
  {
    cout << "\nBoost Tessellator:\n" << endl;
    BoostTessellator<double> tessellator;
    generateMesh(tessellator);
  }
#endif

   // cout << "\nVoro 2D Tessellator:\n" << endl;
   // VoroPP_2d<double> voro;
   // generateMesh(voro);
   // cout << "Voro 2D: PASS" << endl;

  cout << "PASS" << endl;

#if HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
