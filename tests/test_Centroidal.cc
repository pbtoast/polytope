
#include <iostream>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>
#include <sstream>
#include <algorithm>

#include "polytope.hh"
#include "Boundary2D.hh"
#include "Generators.hh"
#include "polytope_test_utilities.hh"
#include "polytope_geometric_utilities.hh"

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
   
   unsigned nPoints = 100;     // Number of generators
   unsigned nIter   = 100;      // Number of iterations

   // Set up boundary and disperse random generator locations
   Boundary2D<double> boundary;
   boundary.setDonut();
   Generators<2,double> generators(boundary);
   generators.randomPoints(nPoints);
   std::vector<double> points;
   for (unsigned i = 0; i != nPoints; ++i) {
     if (boundary.testInside(&generators.mPoints[2*i])) {
       std::copy(&generators.mPoints[2*i], &generators.mPoints[2*i+2], std::back_inserter(points));
     }
   }
   
   // Initialize mesh and tessellator
   Tessellation<2,double> mesh;
   TriangleTessellator<double> triangle;

   for( int iter = 0; iter < nIter; ++iter ){
      mesh.clear();
      triangle.tessellate(points, boundary.mPLCpoints, boundary.mPLC, mesh);
      
#if HAVE_SILO
      vector<double> index(points.size()/2);
      vector<double> genx (points.size()/2);
      vector<double> geny (points.size()/2);
      for (int i = 0; i < mesh.cells.size()/2; ++i){
         index[i] = double(i);
         genx[i]  = points[2*i];
         geny[i]  = points[2*i+1];
      }
      map<string,double*> nodeFields, edgeFields, faceFields, cellFields;
      cellFields["cell_index"   ] = &index[0];
      cellFields["cell_center_x"] = &genx[0];
      cellFields["cell_center_y"] = &geny[0];
      ostringstream os;
      os << "test_Centroidal_" << nPoints << "points_" << iter;
      char dirname[1024];
      snprintf(dirname, 1024, "%s-%d", os.str().c_str(), 0);
      string masterDirName = dirname;
      polytope::SiloWriter<2, double>::write(mesh, nodeFields, edgeFields, 
                                             faceFields, cellFields, os.str());
#endif
      
      for (int i = 0; i < mesh.cells.size(); ++i){
         double tmp[2];
         geometry::computeCellCentroid(mesh, i, tmp);
         points[2*i  ] = 0.5*(points[2*i  ] + tmp[0]);
         points[2*i+1] = 0.5*(points[2*i+1] + tmp[1]);
      }
   }

#if HAVE_MPI
   MPI_Finalize();
#endif
   return 0;
}
