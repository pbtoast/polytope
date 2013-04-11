// test_Circumcenters
//
// Tests for Triangle Tessellator involving degenerate Delaunay
// 1. small-area triangle has a circumcenter approaching infinity

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

using namespace std;
using namespace polytope;

// -----------------------------------------------------------------------
// outputSiloMesh
// -----------------------------------------------------------------------
void outputSiloMesh(Tessellation<2,double>& mesh,
                    std::vector<double>& points,
                    int ntest) {
#if HAVE_SILO
   vector<double> index(mesh.cells.size());
   vector<double> genx (mesh.cells.size());
   vector<double> geny (mesh.cells.size());
   for (int i = 0; i < mesh.cells.size(); ++i){
      index[i] = double(i);
      if(!points.empty()){
         genx[i] = points[2*i];
         geny[i] = points[2*i+1];
      }
   }
   map<string,double*> nodeFields, edgeFields, faceFields, cellFields;
   cellFields["cell_index"   ] = &index[0];
   cellFields["cell_center_x"] = &genx[0];
   cellFields["cell_center_y"] = &geny[0];
   ostringstream os;
   os << "test_Circumcenters";
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

  vector<double> PLCpoints;
  PLCpoints.push_back(0.0);  PLCpoints.push_back(0.0);
  PLCpoints.push_back(3.0);  PLCpoints.push_back(0.0);
  PLCpoints.push_back(3.0);  PLCpoints.push_back(2.0);
  PLCpoints.push_back(0.0);  PLCpoints.push_back(2.0);

  PLC<2,double> boundary;
  boundary.facets.resize(4, vector<int>(2));
  for (unsigned i=0; i != 4; ++i) {
     boundary.facets[i][0] = i;
     boundary.facets[i][1] = (i+1) % 4;
  }
  
  vector<double> points;
  points.push_back(1.5);  points.push_back(1.5);
  points.push_back(0.5);  points.push_back(0.5);
  points.push_back(2.5);  points.push_back(0.5);
  points.push_back(0.5);  points.push_back(1.5);
  points.push_back(1.5);  points.push_back(1.5);
  points.push_back(2.5);  points.push_back(1.5);
  
  TriangleTessellator<double> tessellator;
  //BoostTessellator<double> tessellator;

  int N = 34;
  double vert1[2] = {0.5, 0.5}, vert2[2] = {2.5, 0.5}, circumcenter[2];
  double displacement = 1.0;
  for (int i=0; i != N; ++i) {
     displacement *= 0.5;
     double vert3[2] = {1.5,0.5+displacement};
     geometry::computeCircumcenter2d(vert1, vert2, vert3, circumcenter);
     cout << "Step " << i << ":" << endl;
     cout << "  New generator position  = (1.5," << 0.5+displacement << ")" << endl;
     cout << "  Degenerate tri height   = " << displacement << endl;
     cout << "  Degenerate circumcenter = (" << circumcenter[0] << ","
          << circumcenter[1] << ")" << endl;
     points[1] = 0.5 + displacement;
     Tessellation<2,double> mesh;
     tessellator.tessellate(points,PLCpoints,boundary,mesh);
     //tessellator.tessellate(points,mesh);

     outputSiloMesh(mesh,points,i);

     cout << "          Number of nodes = " << mesh.nodes.size()/2 << endl;
     cout << "          Number of faces = " << mesh.faces.size() << endl << endl;
     // cout << "  Nodes:" << endl;
     // for (int j=0; j != mesh.nodes.size()/2; ++j){
     //    cout << scientific << setprecision(numeric_limits<double>::digits)
     //         <<"    (" << mesh.nodes[2*j] << "," << mesh.nodes[2*j+1] << ")" << endl;
     // }
  }

  cout << "PASS" << endl;

#if HAVE_MPI
   MPI_Finalize();
#endif
   return 0;
}
