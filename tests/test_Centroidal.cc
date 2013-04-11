
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
// lloyd
// -----------------------------------------------------------------------
void lloyd(Tessellation<2,double>& mesh,
           vector<double>& points) {
   for (unsigned i = 0; i < mesh.cells.size(); ++i){
      double cent[2], area;
      geometry::computeCellCentroidAndSignedArea(mesh, i, 1.0e-12, cent, area);
      points[2*i  ] = 0.5*(points[2*i  ] + cent[0]);
      points[2*i+1] = 0.5*(points[2*i+1] + cent[1]);
   }
}

// -----------------------------------------------------------------------
// outputMesh
// -----------------------------------------------------------------------
void outputMesh(Tessellation<2,double>& mesh, 
                const vector<double>& points,
                int nstep) {
#if HAVE_SILO
   vector<double> index(mesh.cells.size());
   vector<double> genx (mesh.cells.size());
   vector<double> geny (mesh.cells.size());
   for (int i = 0; i < mesh.cells.size(); ++i){
      index[i] = double(i);
      genx[i] = points[2*i  ];
      geny[i] = points[2*i+1];
   }
   map<string,double*> nodeFields, edgeFields, faceFields, cellFields;
   cellFields["cell_index"   ] = &index[0];
   cellFields["cell_center_x"] = &genx[0];
   cellFields["cell_center_y"] = &geny[0];
   ostringstream os;
   os << "test_Centroidal";
   polytope::SiloWriter<2, double>::write(mesh, nodeFields, edgeFields, 
                                          faceFields, cellFields, os.str(),
                                          nstep, 0.0);
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
     outputMesh(mesh, points, iter);
     lloyd(mesh,points);
   }
   
   cout << "PASS" << endl;
   
#if HAVE_MPI
   MPI_Finalize();
#endif
   return 0;
}
