// test_TiltedLattice
//
// Initialize generators on an NxN Cartesian lattice for N in [2,100],
// then solidly rotate the lattice by 1 degree.
// Check that the resulting tessellation is still a lattice of quads.
// -----------------------------------------------------------------------
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

#ifdef HAVE_MPI
#include "mpi.h"
#endif

using namespace std;
using namespace polytope;

double tiltX(const double x, const double y, const double xcen, const double ycen, const double theta) {
  return xcen + cos(theta)*(x-xcen) - sin(theta)*(y-ycen);
}

double tiltY(const double x, const double y, const double xcen, const double ycen, const double theta) {
  return ycen + cos(theta)*(y-ycen) + sin(theta)*(x-xcen);
}

// -----------------------------------------------------------------------
// checkCartesianMesh
// -----------------------------------------------------------------------
void checkCartesianMesh(Tessellation<2,double>& mesh, unsigned nx, unsigned ny)
{
  POLY_CHECK(mesh.nodes.size()/2 == (nx + 1)*(ny + 1) );
  POLY_CHECK(mesh.cells.size()   == nx*ny );
  POLY_CHECK(mesh.faces.size()   == nx*(ny + 1) + ny*(nx + 1) );
  for (unsigned i = 0; i != nx*ny; ++i) {
    POLY_CHECK2(mesh.cells[i].size() == 4,
		"Error at cell " << i << ": number of faces = " << mesh.cells[i].size()
		<< endl << endl << mesh);
  }

  vector<set<unsigned> > nodeCells = mesh.computeNodeCells();
  for (unsigned i = 0; i != (nx+1)*(ny+1); ++i) {
    POLY_CHECK2((nodeCells[i].size() == 4) or
		(nodeCells[i].size() == 2) or
		(nodeCells[i].size() == 1),
		"Error at cell " << i << ": number of nodes = " << nodeCells[i].size());
  }
}

// -----------------------------------------------------------------------
// generateMesh
// -----------------------------------------------------------------------
void generateMesh(Tessellator<2,double>& tessellator,
                  const bool checkLattice)
{
  // Constants
  const double degToRad = 2.0*M_PI/360;

  // Parameters
  const double degrees = 4.0;     // Angle of tilt
  const unsigned Nmin = 2;        // Minimum value of N for NxN lattice
  const unsigned Nmax = 100;      // Maximum value of N for NxN lattice
  const unsigned dumpEvery = 10;  // Set to 1 to dump EVERY mesh

  // Derived parameters
  const double angle = degrees * degToRad;

  double boundaryPoints[8] = {0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0};
  vector<double> plcPoints(8);
  PLC<2,double> plc;
  plc.facets.resize(4, vector<int>(2));
  for (unsigned i = 0; i < 4; ++i) {
    plcPoints[2*i  ] = tiltX(boundaryPoints[2*i], boundaryPoints[2*i+1], 0.5, 0.5, angle);
    plcPoints[2*i+1] = tiltY(boundaryPoints[2*i], boundaryPoints[2*i+1], 0.5, 0.5, angle);
    plc.facets[i][0] = i;
    plc.facets[i][1] = (i+1)%4;
  }

  for (unsigned N = Nmin; N < Nmax+1; ++N){
    cout << "Testing N=" << N << endl;

    const double delta = 1.0/N;

    // Create generators
    vector<double> points(2*N*N);
    int index = 0;
    for (unsigned iy = 0; iy < N; ++iy) {
      double yi = (iy + 0.5)*delta;
      for (unsigned ix = 0; ix < N; ++ix) {
	double xi = (ix + 0.5)*delta;
	points[2*index  ] = tiltX(xi, yi, 0.5, 0.5, angle);
	points[2*index+1] = tiltY(xi, yi, 0.5, 0.5, angle);
	index++;
      }
    }

    Tessellation<2,double> mesh;
    tessellator.tessellate(points, plcPoints, plc, mesh);
    if (N % dumpEvery == 0)  outputMesh(mesh, "tiltedLattice_" + tessellator.name(), points, N);
    
    // CHECKS:
    cout << "   num mesh nodes : " << mesh.nodes.size()/2 << endl;
    cout << "   num mesh cells : " << mesh.cells.size()   << endl;
    cout << "   num mesh faces : " << mesh.faces.size()   << endl;
    if (checkLattice)  checkCartesianMesh(mesh, N, N);
  }
}


// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
// -----------------------------------------------------------------------
int main(int argc, char** argv)
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif


#ifdef HAVE_TRIANGLE
  {
    cout << "\nTriangle Tessellator:\n" << endl;
    TriangleTessellator<double> tessellator;
    generateMesh(tessellator, true);
  }
#endif   

#ifdef HAVE_BOOST_VORONOI
  {
    cout << "\nBoost Tessellator:\n" << endl;
    BoostTessellator<double> tessellator;
    generateMesh(tessellator, false);
    
    // NOTE: We do not check to see if BoostTessellator
    // resolves the degenerate lattice mesh. Boost.Voronoi 
    // requires integer input data; it is an unrealistic
    // expectation to think a degenerate Cartesian mesh
    // could be resolved with generators snapped to an
    // integer grid.  -DPS 01/20/2015
  }
#endif

  cout << "PASS" << endl;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
