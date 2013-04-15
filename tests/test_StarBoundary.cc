// test_StarBoundary

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
// outputMesh
// -----------------------------------------------------------------------
void outputMesh(Tessellation<2,double>& mesh, 
                std::vector<double>& points,
                int ntest) {
#if HAVE_SILO
   vector<double> index(mesh.cells.size());
   vector<double> genx (mesh.cells.size());
   vector<double> geny (mesh.cells.size());
   for (int i = 0; i < mesh.cells.size(); ++i){
      index[i] = double(i);
      genx[i]  = points[2*i];
      geny[i]  = points[2*i+1];
   }
   map<string,double*> nodeFields, edgeFields, faceFields, cellFields;
   cellFields["cell_index"   ] = &index[0];
   cellFields["cell_center_x"] = &genx[0];
   cellFields["cell_center_y"] = &geny[0];
   ostringstream os;
   os << "test_StarBoundary_test_" << ntest;
   polytope::SiloWriter<2, double>::write(mesh, nodeFields, edgeFields, 
                                          faceFields, cellFields, os.str());
#endif
}

// -----------------------------------------------------------------------
// test
// -----------------------------------------------------------------------
void test(Tessellator<2,double>& tessellator) {
  unsigned i;
  int test = 1;
  
  // Set the boundary and tessellator
  Boundary2D<double> boundary;
  boundary.setStarWithHole();
  TriangleTessellator<double> triangle;
  
  // 9 input generators points
  double gens[18] = {0.11,  0.11,
                     0.10,  0.44,
                     0.05,  0.70,
                    -0.40,  0.10,
                    -0.48,  0.20,
                    -0.30, -0.32,
                     0.45, -0.55,
                     0.55,  0.05,
                     0.80,  0.20};

  // Test 1: Input generators, no points on boundary
  {
    cout << "\nTest 1: 9 input generators" << endl;
    std::vector<double> points;
    for (i = 0; i < 9; ++i){
      points.push_back(gens[2*i  ]);
      points.push_back(gens[2*i+1]);
    }
    Tessellation<2,double> mesh;
    triangle.tessellate(points,boundary.mPLCpoints,boundary.mPLC,mesh);
    outputMesh(mesh,points,test);
    ++test;
  }
  

  // Test 2: Input generators + boundary generators
  {
    cout << "\nTest 2: 9 input generators + boundary generators" << endl;
    std::vector<double> points;
    for (i = 0; i < 9; ++i){
      points.push_back(gens[2*i  ]);
      points.push_back(gens[2*i+1]);
    }
    points.insert(points.begin(),
                  boundary.mPLCpoints.begin(),
                  boundary.mPLCpoints.end());
    Tessellation<2,double> mesh;
    //triangle.tessellate(points,mesh);
    //triangle.tessellate(points,points,boundary.mPLC,mesh);
    triangle.tessellate(points,boundary.mPLCpoints,boundary.mPLC,mesh);
    outputMesh(mesh,points,test);
    ++test;
  }
  

  // Test 3: 800 random generators
  {
    cout << "\nTest 3: 800 random generators" << endl;
    Generators<2,double> generators(boundary);
    generators.randomPoints(800);
    Tessellation<2,double> mesh;
    triangle.tessellate(generators.mPoints,boundary.mPLCpoints,boundary.mPLC,mesh);
    outputMesh(mesh,generators.mPoints,test);
    ++test;
  }
  

  // Test 4: 2000 random generators
  {
    cout << "\nTest 4: 2000 random generators" << endl;
    Generators<2,double> generators(boundary);
    generators.randomPoints(2000);
    Tessellation<2,double> mesh;
    triangle.tessellate(generators.mPoints,boundary.mPLCpoints,boundary.mPLC,mesh);
    outputMesh(mesh,generators.mPoints,test);
    ++test;
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


#if HAVE_TRIANGLE
  {
    cout << "\nTriangle Tessellator:\n" << endl;
    TriangleTessellator<double> tessellator;
    test(tessellator);  
  }
#endif   


#if HAVE_BOOST_VORONOI
  {
    cout << "\nBoost Tessellator:\n" << endl;
    BoostTessellator<double> tessellator;
    test(tessellator);
  }
#endif
   
  cout << "PASS" << endl;

#if HAVE_MPI
  MPI_Finalize();
#endif
   return 0;
}
