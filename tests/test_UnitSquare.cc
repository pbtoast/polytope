// test_UnitSquare
//
// Mesh a unit square with (nx-by-nx) Cartesian generators for nx in [2,100].
// Perform checks on the resulting tessellation to see if it is indeed Cartesian.
// Both Triangle and Voro++ 2D tessellators are tested here.

#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>
#include <sstream>

#include "polytope.hh"
#include "Boundary2D.hh"
#include "Generators.hh"
#include "polytope_test_utilities.hh"

#include "timingUtilities.hh"

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
   
   const unsigned imin   = 2;
   const unsigned imax   = 101;
   const unsigned iscale = 1;

   const unsigned Imax  = 99;
   const unsigned Nmin  = 2;
   const unsigned scale = 1;
   
   for (unsigned i = 0; i != Imax; ++i) {
     unsigned nx = Nmin + scale*i;
     cout << "Testing nx=" << nx << endl;

      // Create generators
      std::vector<unsigned> nxny(2,nx);
      generators.cartesianPoints(nxny);
      Tessellation<2,double> mesh;
      
//      ifstream os;
//      double x, y, z;
//      vector<double> points;
//      os.open("oneMillionRandomGenerators.txt");
//      for (unsigned ii = 0; ii != 1000000; ++ii) {
//        os >> x >> y >> z;
//        points.push_back(x);
//        points.push_back(y);
//      }
//      os.close();

     Timing::Time start = Timing::currentTime();
     tessellator.tessellate(generators.mPoints, boundary.mPLCpoints, boundary.mPLC, mesh);
     //tessellator.tessellate(points, boundary.mPLCpoints, boundary.mPLC, mesh);
     double time = Timing::difference(start, Timing::currentTime());
     
     cout << "Mesh construction took " << time << " seconds." << endl;
     
     // CHECKS:
     cout << "   num mesh nodes : " << mesh.nodes.size()/2 << endl;
     cout << "   num mesh cells : " << mesh.cells.size()   << endl;
     cout << "   num mesh faces : " << mesh.faces.size()   << endl;
     checkCartesianMesh(mesh,nx,nx);

//      ofstream os;
//      os.open("oneMillionRandomGenerators.txt");
//      os << "# x \t\t y \t\t z" << endl;
//      for (unsigned j = 0; j != generators.mPoints.size()/2; ++j) {
//        os << generators.mPoints[2*j] << " \t " << generators.mPoints[2*j+1] << " \t " << 0.0 << endl;
//      }
//      os.close();

//      unsigned ic = 692152;
//      cout << "Cell " << ic << " has face indices" << endl;
//      for (vector<int>::iterator itr = mesh.cells[ic].begin();
// 	  itr != mesh.cells[ic].end(); ++itr) {
//        unsigned iface = (*itr < 0) ? ~(*itr) : *itr;
//        cout << "   " << iface << " with nodes" << endl;
//        unsigned inode0 = (*itr > 0) ? mesh.faces[iface][0] : mesh.faces[iface][1];
//        unsigned inode1 = (*itr > 0) ? mesh.faces[iface][1] : mesh.faces[iface][0];
//        cout << "      (" << mesh.nodes[2*inode0] << "," << mesh.nodes[2*inode0+1] << ")" << endl
// 	    << "      (" << mesh.nodes[2*inode1] << "," << mesh.nodes[2*inode1+1] << ")" << endl;
//        unsigned c0 = mesh.faceCells[iface][0] < 0 ? ~mesh.faceCells[iface][0] : mesh.faceCells[iface][0];
//        unsigned c1 = mesh.faceCells[iface][1] < 0 ? ~mesh.faceCells[iface][1] : mesh.faceCells[iface][1];
//        unsigned icell = (c0 == ic) ? c1 : c0;
//        cout << "   and neighboring cell " << icell << endl;
//      }
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
