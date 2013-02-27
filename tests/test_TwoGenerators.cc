// test_twoGenerators
//
// If coordMax in TriangleTessellator is set too low, then a wrap-around effect
// will become visible when trying to tessellate a unit-square boundary

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
// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
int main(int argc, char** argv)
{
#if HAVE_MPI
   MPI_Init(&argc, &argv);
#endif

  Boundary2D<double> boundary;
  boundary.unitSquare();
  Generators<2,double> generators(boundary);
  
  double point1[2] = {-0.25, -0.125};  generators.addGenerator(point1);
  double point2[2] = { 0.25,  0.125};  generators.addGenerator(point2);
  
  Tessellation<2,double> mesh;
  TriangleTessellator<double> triangle;
  triangle.tessellate( generators.mPoints, boundary.mPLCpoints, boundary.mPLC, mesh );
  
  cerr << "\n\n\n\n\n" << mesh << endl;

#if HAVE_SILO
  vector<double> index( mesh.cells.size());
  for (int i = 0; i < mesh.cells.size(); ++i) index[i] = double(i);
  map<string,double*> nodeFields, edgeFields, faceFields, cellFields;
  cellFields["cell_index"] = &index[0];
  ostringstream os;
  os << "test_twoGenerators_mesh";
  polytope::SiloWriter<2, double>::write(mesh, nodeFields, edgeFields, 
                                         faceFields, cellFields, os.str());
#endif
  
  POLY_CHECK(mesh.nodes.size()/2 == 6);
  POLY_CHECK(mesh.cells.size()   == 2);
  POLY_CHECK(mesh.faces.size()   == 7);
  for (unsigned i = 0; i != 2; ++i) POLY_CHECK(mesh.cells[i].size() == 4);
  const double area = computeTessellationArea( mesh );
  const double fracerr = std::abs(boundary.mArea - area)/boundary.mArea;
  const double tol = 1.0e-6;
  POLY_CHECK2( fracerr < tol, ""
               << "              Area  = " << area << endl
               << "              Error = " << boundary.mArea - area << endl
               << "   Fractional error = " << fracerr << endl );

  cout << "PASS" << endl;

#if HAVE_MPI
   MPI_Finalize();
#endif
   return 0;
}
