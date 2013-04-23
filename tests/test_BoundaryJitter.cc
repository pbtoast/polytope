// test_BoundaryJitter
//
// Set up a 10x10 set of Cartesian generators. Hold the 4 corners fixed
// and the 8x8 interior generators fixed. Repeatedly jitter the remaining
// generators on the boundary by random 10^-16 perturbations and see if
// near-infinite circumcenters cause us to fail

#include <iostream>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>
#include <sstream>

#include "polytope.hh"
#include "Boundary2D.hh"
#include "polytope_test_utilities.hh"

#if HAVE_MPI
#include "mpi.h"
#endif

#define POLY_CHECK_BOOL(x) if (!(x)) { return false; }

using namespace std;
using namespace polytope;

// -----------------------------------------------------------------------
// minLength
// -----------------------------------------------------------------------
double minLength(Tessellation<2,double>& mesh)
{
   double faceLength = FLT_MAX;
   for (unsigned iface = 0; iface != mesh.faces.size(); ++iface)
   {
      POLY_ASSERT( mesh.faces[iface].size() == 2 );
      const unsigned inode0 = mesh.faces[iface][0];
      const unsigned inode1 = mesh.faces[iface][1];
      double x0 = mesh.nodes[2*inode0], y0 = mesh.nodes[2*inode0+1];
      double x1 = mesh.nodes[2*inode1], y1 = mesh.nodes[2*inode1+1];
      double len = (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0);
      faceLength = min( faceLength, sqrt(len) );
   }
   return faceLength;
}

// -----------------------------------------------------------------------
// checkIfCartesian
// -----------------------------------------------------------------------
bool checkIfCartesian(Tessellation<2,double>& mesh, unsigned nx, unsigned ny) {
   POLY_CHECK_BOOL(mesh.nodes.size()/2 == (nx + 1)*(ny + 1) );
   POLY_CHECK_BOOL(mesh.cells.size()   == nx*ny );
   POLY_CHECK_BOOL(mesh.faces.size()   == nx*(ny + 1) + ny*(nx + 1) );
   for (unsigned i = 0; i != nx*ny; ++i) POLY_CHECK_BOOL(mesh.cells[i].size() == 4);
   
   std::vector<std::set<unsigned> > nodeCells = mesh.computeNodeCells();
   for (unsigned i = 0; i != (nx+1)*(ny+1); ++i) {
      POLY_CHECK_BOOL( (nodeCells[i].size() == 4) ||
                       (nodeCells[i].size() == 2) ||
                       (nodeCells[i].size() == 1) );
   }
   return true;
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
   os << "test_BoundaryJitter";
   polytope::SiloWriter<2, double>::write(mesh, nodeFields, edgeFields, 
                                          faceFields, cellFields, os.str(),
                                          nstep, 0.0);
#endif
}


// -----------------------------------------------------------------------
// computeJitterMask
// -----------------------------------------------------------------------
vector<unsigned> computeJitterMask(const unsigned nx) {
   const unsigned numGenerators = nx*nx;
   vector<unsigned> jitterMask(numGenerators, 0);
   unsigned i1, i2, i3, i4;
   for (unsigned i = 1; i != nx-1; ++i) {
      i1 = i;
      i2 = i*nx;
      i3 = (i+1)*nx;
      i4 = nx*(nx-1)+i;
      POLY_ASSERT(i1 < numGenerators and i2 < numGenerators and
                  i3 < numGenerators and i4 < numGenerators);
      jitterMask[i1] = 1;
      jitterMask[i2] = 1;
      jitterMask[i3] = 1;
      jitterMask[i4] = 1;
   }
   return jitterMask;
}
   

// -----------------------------------------------------------------------
// jitterPoints
// -----------------------------------------------------------------------
void jitterPoints(vector<double>& points,
                  vector<unsigned>& jitterMask,
                  const double epsilon) {
  const unsigned numGenerators = points.size()/2;
  POLY_ASSERT(numGenerators == jitterMask.size());
  for (unsigned i = 0; i != numGenerators; ++i) {
    if (jitterMask[i] == 1) {
      points[2*i  ] += epsilon*(random01()-0.5);
      points[2*i+1] += epsilon*(random01()-0.5);
    }
  }
}


// -----------------------------------------------------------------------
// test
// -----------------------------------------------------------------------
void test(Tessellator<2,double>& tessellator) {
  // Input generator parameters
  const unsigned nx = 20;
  const double xmin = -0.5, xmax = 0.5;
  const double ymin = -0.5, ymax = 0.5;
  const double dx = (xmax-xmin)/nx,  dy = (ymax-ymin)/nx;
  POLY_ASSERT(nx > 2);

  // Jitter factor
  const double epsilon = 2.0e-10;

  // Set the boundary
  Boundary2D<double> boundary;
  boundary.setUnitSquare();
  vector<unsigned> jitterMask = computeJitterMask(nx);

  // Initialize the mesh
  Tessellation<2, double> mesh;

  // The initial generators
  vector<double> points;
  unsigned ix, iy;
  double xi, yi;
  for (iy = 0; iy != nx; ++iy) {
    yi = ymin + (iy + 0.5)*dy;
    for (ix = 0; ix != nx; ++ix) {
      xi = xmin + (ix + 0.5)*dx;
      points.push_back(xi);
      points.push_back(yi);
    }
  }
    
  for (unsigned i = 0; i != 100; ++i) {
    mesh.clear();
    jitterPoints(points, jitterMask, epsilon);
    tessellator.tessellate(points, boundary.mPLCpoints, boundary.mPLC, mesh);
    outputMesh(mesh, points, i);
    bool isCartesian = checkIfCartesian(mesh,nx,nx);
    if(!isCartesian) 
       cout << "Degeneracy threshold reached! Minimum face length = " << minLength(mesh) << endl;
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
