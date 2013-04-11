// test_MeshEditor
//
// Unit tests for the MeshEditor class and cleanEdges routine

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
// outputSiloMesh
// -----------------------------------------------------------------------
void outputSiloMesh(Tessellation<2,double>& mesh,
                    int ntest) {
#if HAVE_SILO
   vector<double> index(mesh.cells.size());
   for (int i = 0; i < mesh.cells.size(); ++i){
      index[i] = double(i);
   }
   map<string,double*> nodeFields, edgeFields, faceFields, cellFields;
   cellFields["cell_index"   ] = &index[0];
   ostringstream os;
   os << "test_MeshEditor";
   polytope::SiloWriter<2, double>::write(mesh, nodeFields, edgeFields, 
                                          faceFields, cellFields, os.str(),
                                          ntest, 0.0);
#endif
}



// -----------------------------------------------------------------------
// main
// -----------------------------------------------------------------------
int main(int argc, char** argv)
{
#if HAVE_MPI
   MPI_Init(&argc, &argv);
#endif

   int itest = 1;

   // Test 1: Single cell, one small face
   {
      const double edgeTol = 0.001;
      cout << "\nTest 1: Eliminate small face on a single 2D cell" << endl;
      Tessellation<2,double> mesh;
      const unsigned ncells = 1;
      const unsigned nfaces = 6;
      const unsigned nnodes = 6;
      mesh.cells.resize(ncells, std::vector<int>(nfaces));
      mesh.faces.resize(nfaces, std::vector<unsigned>(2));
      mesh.faceCells.resize(nfaces, std::vector<int>(1));
      for (int i = 0; i != nfaces; ++i) {
         mesh.cells[0][i] = i;
         mesh.faces[i][0] = i;
         mesh.faces[i][1] = (i+1) % nnodes;
         mesh.faceCells[i][0] = 0;
      }
      double nodePoints[12] = {0.00, 0.0,
                               0.75, 0.0,
                               1.00, 0.5,
                               1.00, 0.50000001,
                               0.75, 1.0,
                               0.00, 1.0};
      mesh.nodes.resize(2*nnodes);
      for (int i = 0; i != 2*nnodes; ++i) mesh.nodes[i] = nodePoints[i];
      MeshEditor<2,double> meshEditor(mesh);
      cout << mesh << endl;
      outputSiloMesh(mesh, itest);
      meshEditor.cleanEdges(edgeTol);
      cout << mesh << endl;
      outputSiloMesh(mesh, itest);
      ++itest;
   }

   cout << "PASS" << endl;

#if HAVE_MPI
   MPI_Finalize();
#endif
   return 0;
}
