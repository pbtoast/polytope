// test_OrphanCases
//
// A collection of difficult test cases for the orphaned cell algorithm.
// All tests involve meshing a circular region with a star-shaped hole
// in the middle with only 20 points. Seeding the random number generator
// provides the input generator locations. Test cases were found through 
// trial-and-error (though with obnoxious frequency!).

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
   os << "test_OrphanCases_test_" << ntest;
   polytope::SiloWriter<2, double>::write(mesh, nodeFields, edgeFields, 
                                          faceFields, cellFields, os.str());
#endif
}

// -----------------------------------------------------------------------
// outputMesh
// -----------------------------------------------------------------------
void outputMesh(Boundary2D<double>& boundary,
                std::vector<double>& points,
                Tessellation<2,double>& mesh, 
                int ntest) {
   POLY_ASSERT( mesh.cells.size() == 20 );
   double area = computeTessellationArea(mesh);
   cout << "Tessellation Area = " << area << endl;
   cout << "Relative error    = " << (boundary.mArea-area)/boundary.mArea << endl;
   outputSiloMesh(mesh, points, ntest);
}

// -----------------------------------------------------------------------
// main
// -----------------------------------------------------------------------
int main(int argc, char** argv)
{
#if HAVE_MPI
   MPI_Init(&argc, &argv);
#endif

   // Initialize boundary and tessellator
   Boundary2D<double> boundary;
   TriangleTessellator<double> triangle;
   
   // Circular region with star-shaped hole
   int bType = 5;
   boundary.setDefaultBoundary(bType);

   int i = 1;

   // Test 1: Cell parents multiple orphans
   {
      srand(10489593);
      cout << "\nTest 1: Cell parents multiple orphans" << endl;
      Generators<2,double> generators(boundary);
      generators.randomPoints(20);
      Tessellation<2,double> mesh;
      tessellate2D(generators.mPoints, boundary, triangle, mesh);
      outputMesh(boundary, generators.mPoints, mesh, i);
      ++i;
   }

   // Test 2: Orphan neighbors are also parents of orphans
   {
      srand(10489594);
      cout << "\nTest 2: Orphan neighbors are also parents of orphans" << endl;
      Generators<2,double> generators(boundary);
      generators.randomPoints(20);
      Tessellation<2,double> mesh;
      tessellate2D(generators.mPoints, boundary, triangle, mesh);
      outputMesh(boundary, generators.mPoints, mesh, i);
      ++i;
   }

   // Test 3: Overlapping orphans
   {
      srand(10489609);
      cout << "\nTest 3: Overlapping orphans" << endl;
      Generators<2,double> generators(boundary);
      generators.randomPoints(20);
      Tessellation<2,double> mesh;
      tessellate2D(generators.mPoints, boundary, triangle, mesh);
      outputMesh(boundary, generators.mPoints, mesh, i);
      ++i;
   }

   // Test 4: Empty orphan neighbor set
   {
      srand(10489611);
      cout << "\nTest 4: Empty orphan neighbor set" << endl;
      Generators<2,double> generators(boundary);
      generators.randomPoints(20);
      Tessellation<2,double> mesh;
      tessellate2D(generators.mPoints, boundary, triangle, mesh);
      outputMesh(boundary, generators.mPoints, mesh, i);
      ++i;
   }

   // Test 5: Boost.Geometry calls invalid overlay exception
   {
      srand(10489612);
      cout << "\nTest 5: Boost.Geometry calls invalid overlay exception" << endl;
      Generators<2,double> generators(boundary);
      generators.randomPoints(20);
      Tessellation<2,double> mesh;
      tessellate2D(generators.mPoints, boundary, triangle, mesh);
      outputMesh(boundary, generators.mPoints, mesh, i);
      ++i;
   }

   // Test 6: Nonconvex boundary with three internal generators
   {
      cout << "\nTest 6: Nonconvex boundary with three internal generators" << endl;
      std::vector<double> points;
      points.push_back(0.05); points.push_back(0.60);
      points.push_back(0.40); points.push_back(0.10);
      points.push_back(0.60); points.push_back(0.50);
      std::vector<double> PLCpoints;
      PLCpoints.push_back(0.0); PLCpoints.push_back(0.0);
      PLCpoints.push_back(0.1); PLCpoints.push_back(0.0);
      PLCpoints.push_back(0.2); PLCpoints.push_back(0.8);
      PLCpoints.push_back(0.3); PLCpoints.push_back(0.0);
      PLCpoints.push_back(1.0); PLCpoints.push_back(0.0);
      PLCpoints.push_back(1.0); PLCpoints.push_back(1.0);
      PLCpoints.push_back(0.0); PLCpoints.push_back(1.0);
      PLC<2,double> boundary;
      boundary.facets.resize(7, std::vector<int>(2));
      for (int j = 0; j < 7; ++j){
         boundary.facets[j][0] = j;
         boundary.facets[j][1] = (j+1) % 7;
      }
      Tessellation<2,double> mesh;
      triangle.tessellate(points, PLCpoints, boundary, mesh);
      outputSiloMesh(mesh, points, i);
      ++i;
   }

   cout << "PASS" << endl;

#if HAVE_MPI
   MPI_Finalize();
#endif
   return 0;
}
