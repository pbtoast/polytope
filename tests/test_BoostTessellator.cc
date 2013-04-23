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
// test1
//------------------------------------------------------------------------
void test1(Tessellator<2,double>& tessellator) {
  
  // Output name
  string testName = "BoostTessellator_test1";

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

  // Tessellate unbounded
  tessellator.tessellate(points, mesh);
  outputMesh(mesh,testName,points,1);
  mesh.clear();

  // Tessellate bounded
  tessellator.tessellate(points, PLCpoints, boundary, mesh);
  outputMesh(mesh,testName,points,2);
}


//------------------------------------------------------------------------
// test2
//------------------------------------------------------------------------
void test2(Tessellator<2,double>& tessellator) {
  
  // Output name
  string testName = "BoostTessellator_test2";

  // Rectangular boundary
  vector<double> PLCpoints;
  PLCpoints.push_back(0.0);  PLCpoints.push_back(0.0);
  PLCpoints.push_back(1.0);  PLCpoints.push_back(0.0);
  PLCpoints.push_back(1.0);  PLCpoints.push_back(0.5);
  PLCpoints.push_back(0.0);  PLCpoints.push_back(0.5);

  // Three generators
  vector<double> points;
  points.push_back(0.2);  points.push_back(0.25   );
  points.push_back(0.6);  points.push_back(0.25   );
  points.push_back(0.8);  points.push_back(0.25001);
  //points.push_back(0.8);  points.push_back(0.24999);

  // Facets
  PLC<2,double> boundary;
  int nSides = 4;
  boundary.facets.resize(nSides, vector<int>(2));
  for (int i = 0; i != nSides; ++i) {
    boundary.facets[i][0] = i;
    boundary.facets[i][1] = (i+1)%nSides;
  }
  
  Tessellation<2,double> mesh;

  // Tessellate unbounded
  tessellator.tessellate(points, mesh);
  outputMesh(mesh,testName,points,1);
  mesh.clear();

  // Tessellate bounded
  tessellator.tessellate(points, PLCpoints, boundary, mesh);
  outputMesh(mesh,testName,points,2);

  //cout << mesh << endl;
}


//------------------------------------------------------------------------
// main
//------------------------------------------------------------------------
int 
main(int argc, char** argv) 
{
#if HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  BoostTessellator<double> tessellator;
  
  {
    cout << "\nTest 1" << endl;
    test1(tessellator);
  }

  {
    cout << "\nTest 2" << endl;
    test2(tessellator);
  }


  cout << "PASS" << endl;

#if HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
//------------------------------------------------------------------------
