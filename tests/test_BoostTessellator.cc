#include <iostream>
#include <cmath>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>

#include "polytope.hh"
#include "polytope_test_utilities.hh"

#if HAVE_MPI
#include "mpi.h"
#endif

using namespace std;
using namespace polytope;

//------------------------------------------------------------------------
// main
//------------------------------------------------------------------------
int 
main(int argc, char** argv) 
{
#if HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  // Output name
  string testName = "BoostTessellator";

  // Create the generators.
  vector<double> points;
  const int nx = 3;
  const double x1 = 0.0, y1 = 0.0;
  const double x2 = 1.0, y2 = 1.0;
  const double dx = (x2 - x1)/nx, dy = (y2 - y1)/nx;
  double low [2] = { numeric_limits<double>::max(),  numeric_limits<double>::max()};
  double high[2] = {-numeric_limits<double>::max(), -numeric_limits<double>::max()};
  unsigned ix, iy;
  double xi, yi;
  for (iy = 0; iy != nx; ++iy) {
    yi = y1 + (iy + 0.5)*dy;
    for (ix = 0; ix != nx; ++ix) {
      xi = x1 + (ix + 0.5)*dx;
      points.push_back(xi);
      points.push_back(yi);
      low [0] = min(low [0], xi);
      low [1] = min(low [1], yi);
      high[0] = max(high[0], xi);
      high[1] = max(high[1], yi);
    }
  }

  vector<double> PLCpoints;
  PLCpoints.push_back(x1);  PLCpoints.push_back(y1);
  PLCpoints.push_back(x2);  PLCpoints.push_back(y1);
  PLCpoints.push_back(x2);  PLCpoints.push_back(y2);
  PLCpoints.push_back(x1);  PLCpoints.push_back(y2);

  PLC<2,double> boundary;
  int nSides = 4;
  boundary.facets.resize(nSides, vector<int>(2));
  for (int i = 0; i != nSides; ++i) {
    boundary.facets[i][0] = i;
    boundary.facets[i][1] = (i+1)%nSides;
  }

  Tessellation<2,double> mesh;
  BoostTessellator<double> tessellator;
  //tessellator.tessellate(points, mesh);
  tessellator.tessellate(points, PLCpoints, boundary, mesh);
  cout << mesh << endl;

  outputMesh(mesh,testName,points);


  cout << "PASS" << endl;

#if HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
//------------------------------------------------------------------------
