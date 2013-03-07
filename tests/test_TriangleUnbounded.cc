// 

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
// checkMesh
// -----------------------------------------------------------------------
void checkMesh(const Tessellation<2,double>& mesh,
               const unsigned ncells,
               const unsigned nnodes,
               const unsigned nfaces,
               const unsigned ninfNodes,
               const unsigned ninfFaces) {
  POLY_CHECK(mesh.cells.size()   == ncells);
  POLY_CHECK(mesh.nodes.size()/2 == nnodes);
  POLY_CHECK(mesh.faces.size()   == nfaces);
  unsigned infNodeCount=0;
  for(unsigned i = 0; i < mesh.infNodes.size(); ++i){
    if( mesh.infNodes[i] == 1 ) ++infNodeCount;
  }
  POLY_CHECK(infNodeCount == ninfNodes);
  unsigned infFaceCount=0;
  for(unsigned i = 0; i < mesh.infFaces.size(); ++i){
    if( mesh.infFaces[i] == 1 ) ++infFaceCount;
  }
  POLY_CHECK(infFaceCount == ninfFaces);
}

// -----------------------------------------------------------------------
// outputMesh
// -----------------------------------------------------------------------
void outputMesh(Tessellation<2,double>& mesh, int ntest) {
#if HAVE_SILO
   vector<double> index(mesh.cells.size());
   vector<double> genx (mesh.cells.size());
   vector<double> geny (mesh.cells.size());
   for (int i = 0; i < mesh.cells.size(); ++i){
      index[i] = double(i);
   }
   map<string,double*> nodeFields, edgeFields, faceFields, cellFields;
   cellFields["cell_index"   ] = &index[0];
   ostringstream os;
   os << "test_TriangleUnbounded_test_" << ntest;
   polytope::SiloWriter<2, double>::write(mesh, nodeFields, edgeFields, 
                                          faceFields, cellFields, os.str());
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

  int test = 1;
  TriangleTessellator<double> triangle;

  Tessellation<2,double> mesh;

  // Test 1: Circle of generators
  {
    cout << "\nTest 1: Circle of generators" << endl;
    int N = 18;
    vector<double> points(2*N);
    for (unsigned i = 0; i < N; ++i){
       double theta = 2.0*M_PI*i/(N+1);
       points[2*i] = cos(theta);  points[2*i+1] = sin(theta);
    }
    Tessellation<2,double> mesh;
    triangle.tessellate(points, mesh);
    outputMesh(mesh, test);
    checkMesh(mesh,N,N+1,2*N,N,N);
    ++test;
  }

  // Test 2: Two uniform rows of generators
  {
    cout << "\nTest 2: Two uniform rows of generators" << endl;
    int N = 10;
    vector<double> points(4*N);
    for (unsigned i = 0; i < N; ++i){
       points[4*i  ] = i;  points[4*i+1] = -1.0;
       points[4*i+2] = i;  points[4*i+3] =  1.0;
    }
    Tessellation<2,double> mesh;
    triangle.tessellate(points, mesh);
    outputMesh(mesh, test);
    checkMesh(mesh, 2*N, 3*N-1, 5*N-2, 2*N, 2*N);
    ++test;
  }

  // Test 3: Collinear generators with one non-collinear
  {
    cout << "\nTest 3: Collinear generators, except one" << endl;
    int N = 10;
    vector<double> points(2*N);
    for (unsigned i = 0; i < N; ++i)  points[2*i] = double(i);
    points.push_back(4.5);
    points.push_back(1.0);
    Tessellation<2,double> mesh;
    triangle.tessellate(points, mesh);
    outputMesh(mesh, test);
    ++test;
  }

  // Test 4: 2x2 Cartesian Generators
  {
    cout << "\nTest 4: 2x2 Cartesian generators" << endl;
    vector<double> points;
    int N = 2;
    points.push_back(0.0); points.push_back(0.0);
    points.push_back(1.0); points.push_back(0.0);
    points.push_back(1.0); points.push_back(1.0);
    points.push_back(0.0); points.push_back(1.0);
    Tessellation<2,double> mesh;
    triangle.tessellate(points, mesh);
    outputMesh(mesh, test);
    checkMesh(mesh, 4, 5, 8, 4, 4);
    ++test;
  }

  // Test 5: Two generators
  {
    cout << "\nTest 5: Two generators" << endl;
    vector<double> points;
    points.push_back(0.0); points.push_back(0.0);
    points.push_back(1.0); points.push_back(0.0);
    Tessellation<2,double> mesh;
    triangle.tessellate(points, mesh);
    outputMesh(mesh, test);
    checkMesh(mesh, 2, 4, 5, 4, 4);
    ++test;
  }

  // Test 6: Line of generators, uniform
  {
    cout << "\nTest 6: Uniform line of generators" << endl;
    int N=10;
    vector<double> points(2*N, 0);
    for (unsigned i = 0; i < N; ++i)  points[2*i] = double(i);
    Tessellation<2,double> mesh;
    triangle.tessellate(points, mesh);
    outputMesh(mesh, test);
    checkMesh(mesh, N, 2*N, 3*N-1, 2*N, 2*N);
    ++test;
  }

  // Test 7: Line of generators, non-uniform
  {
    cout << "\nTest 7: Non-uniform line of generators" << endl;
    int N = 10;
    vector<double> points(2*N);
    for (unsigned i = 0; i < N; ++i)  points[2*i] = double(i) + random01() - 0.5;
    Tessellation<2,double> mesh;
    triangle.tessellate(points, mesh);
    outputMesh(mesh, test);
    checkMesh(mesh, N, 2*N, 3*N-1, 2*N, 2*N);
    ++test;
  }

  // Test 8: Line of generators, non-uniform, shuffled
  {
    cout << "\nTest 8: Non-uniform line of generators, shuffled" << endl;
    int N = 10;
    vector<double> points(2*N);
    int indices[10] = {5,7,1,2,8,0,3,9,4,6};
    for (unsigned i = 0; i < N; ++i)  points[2*i] = double(indices[i]) + random01() - 0.5;
    Tessellation<2,double> mesh;
    triangle.tessellate(points, mesh);
    outputMesh(mesh, test);
    checkMesh(mesh, N, 2*N, 3*N-1, 2*N, 2*N);    
    ++test;
  }


  cout << "PASS" << endl;

#if HAVE_MPI
   MPI_Finalize();
#endif
   return 0;
}
