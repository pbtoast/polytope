// 

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

  unsigned nPoints = 5;
  double pts[10] = {0.0, 0.0,
                    1.0, 0.0,
                    1.0, 1.0,
                    0.0, 1.0,
                    0.2, 0.5};

  std::vector<double> points;
  for (unsigned i=0; i<nPoints; ++i){
     points.push_back(pts[2*i  ]);
     points.push_back(pts[2*i+1]);
  }
  
  Tessellation<2,double> mesh;
  TriangleTessellator<double> triangle;
  triangle.tessellate( points, mesh );
  
#if HAVE_SILO
  vector<double> index( mesh.cells.size());
  for (int i = 0; i < mesh.cells.size(); ++i) index[i] = double(i);
  map<string,double*> nodeFields, edgeFields, faceFields, cellFields;
  cellFields["cell_index"] = &index[0];
  ostringstream os;
  os << "test_unbounded_mesh";
  polytope::SiloWriter<2, double>::write(mesh, nodeFields, edgeFields, 
                                         faceFields, cellFields, os.str());
#endif
  
  cout << "PASS" << endl;

#if HAVE_MPI
   MPI_Finalize();
#endif
   return 0;
}
