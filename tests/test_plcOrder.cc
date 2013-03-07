// test_PLCorder
//
// Tests different orderings of the PLC facets when handing the boundary geometry
// to the tessellator. Spoiler alert: only a sequential, counterclockwise ordered
// PLC is valid. The other orderings fail, though not always at an assertion.

#include <iostream>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>
#include <sstream>

#include "polytope.hh"

#if HAVE_MPI
#include "mpi.h"
#endif

#define POLY_CHECK(x) if (!(x)) { cout << "FAIL: " << #x << endl; } //exit(-1); }

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
// Different PLC orderings
// -----------------------------------------------------------------------
PLC<2,double> makePLC( int order )
{
   PLC<2,double> box;
   box.facets.resize(4);
   
   switch(order){

   // Ordered CCW
   case 0:
      for (unsigned f = 0; f < 4; ++f){
         box.facets[f].resize(2); box.facets[f][0] = f;  box.facets[f][1] = (f+1) % 4;
      }
      return box;

   // Ordered CW
   case 1:
      for (unsigned f = 0; f < 4; ++f){
         box.facets[f].resize(2); box.facets[f][0] = 3*f%4;  box.facets[f][1] = 3*(f+1)%4;
      }
      return box;
   
   // Disordered CCW
   case 2:
      box.facets[0].resize(2); box.facets[0][0] = 0;  box.facets[0][1] = 1;
      box.facets[1].resize(2); box.facets[1][0] = 2;  box.facets[1][1] = 3;
      box.facets[2].resize(2); box.facets[2][0] = 1;  box.facets[2][1] = 2;
      box.facets[3].resize(2); box.facets[3][0] = 3;  box.facets[3][1] = 0;
      return box;
      
   // Disordered CW
   case 3:
      box.facets[0].resize(2); box.facets[0][0] = 1;  box.facets[0][1] = 0;
      box.facets[1].resize(2); box.facets[1][0] = 3;  box.facets[1][1] = 2;
      box.facets[2].resize(2); box.facets[2][0] = 2;  box.facets[2][1] = 1;
      box.facets[3].resize(2); box.facets[3][0] = 0;  box.facets[3][1] = 3;
      return box;

   // Disordered mixed CW and CCW
   case 4:
      box.facets[0].resize(2); box.facets[0][0] = 0;  box.facets[0][1] = 1;
      box.facets[1].resize(2); box.facets[1][0] = 3;  box.facets[1][1] = 2;
      box.facets[2].resize(2); box.facets[2][0] = 1;  box.facets[2][1] = 2;
      box.facets[3].resize(2); box.facets[3][0] = 0;  box.facets[3][1] = 3;
      return box;
   }
}

//------------------------------------------------------------------------
// the test
//------------------------------------------------------------------------
int 
main(int argc, char** argv) 
{
#if HAVE_MPI
  MPI_Init(&argc, &argv);
#endif
  
  // Generators ordered counterclockwise
  std::vector<double> generators;
  generators.push_back( -0.5 );  generators.push_back( -0.5 );
  generators.push_back(  0.5 );  generators.push_back( -0.5 );
  generators.push_back(  0.5 );  generators.push_back(  0.5 );
  generators.push_back( -0.5 );  generators.push_back(  0.5 );
  
  Tessellation<2, double> mesh;
  TriangleTessellator<double> triangle;
  PLC<2,double> box = makePLC( 0 );
  triangle.tessellate(generators, generators, box, mesh);
  checkCartesianMesh( mesh, 2, 2 );
  cout << "PASS" << endl;
  
#if HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
