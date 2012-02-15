// Try tessellating a simple lattice of generators in a box in parallel.
// We use randomly chosen seed locations to divide up the generators
// between processors.

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <limits>
#include <sstream>
#include <ctime>

#include "polytope.hh"
#include "convexHull_3d.hh"

#define CHECK(x) if (!(x)) { cout << "FAIL: " << #x << endl; exit(-1); }

using namespace std;

//------------------------------------------------------------------------------
// Define our own local random number generator wrapping the standard srand &
// rand methods.
//------------------------------------------------------------------------------
double random01() {
  return double(rand())/RAND_MAX;
}

//------------------------------------------------------------------------------
// Check if a point (p) is contained in a convex PLC surface.
// Returns the facet id for the first facet the point is outside, or the number
// of facets if it's contained.
//------------------------------------------------------------------------------
template<typename RealType>
unsigned
convexContains(const polytope::PLC<3, RealType>& surface,
               const RealType* points,
               const RealType* p,
               const RealType tolerance) {
  bool containmentTest = true;
  unsigned ifacet = 0, n, iv, i, j, k;
  RealType xn, yn, zn, dx1, dy1, dz1, dx2, dy2, dz2, norm;
  while (containmentTest and ifacet < surface.facets.size()) {
    n = surface.facets[ifacet].size();
    ASSERT(n >= 3);
    xn = 0.0; yn = 0.0; zn = 0.0;

    // Compute the facet normal.
    i = surface.facets[ifacet][0];
    for (iv = 1; iv != n - 1; ++iv) {
      j = surface.facets[ifacet][iv];
      k = surface.facets[ifacet][iv + 1];
      dx1 = points[3*j]     - points[3*i];
      dy1 = points[3*j + 1] - points[3*i + 1];
      dz1 = points[3*j + 2] - points[3*i + 2];
      dx2 = points[3*k]     - points[3*i];
      dy2 = points[3*k + 1] - points[3*i + 1];
      dz2 = points[3*k + 2] - points[3*i + 2];
      xn += dy1*dz2 - dz1*dy2;
      yn += dz1*dx2 - dx1*dz2;
      zn += dx1*dy2 - dy1*dx2;
    }
    norm = sqrt(xn*xn + yn*yn + zn*zn);
    xn /= norm;
    yn /= norm;
    zn /= norm;

    // Now test the containment.
    dx1 = p[0] - points[3*i];
    dy1 = p[1] - points[3*i + 1];
    dz1 = p[2] - points[3*i + 2];
    containmentTest = ((dx1*xn + dy1*yn + dz1*zn) < tolerance);
    if (containmentTest) {
      ++ifacet;
    } else {
      cerr << "Facet distance : " << dx1*xn + dy1*yn + dz1*zn << endl
           << "Facet point    : " << points[3*i] << " " << points[3*i + 1] << " " << points[3*i + 2] << endl
           << "Facet normal   : " << xn << " " << yn << " " << zn << endl
           << "Test point     : " << p[0] << " " << p[1] << " " << p[2] << endl;
    }
  }

  return ifacet;
}

//------------------------------------------------------------------------------
// The test itself.
//------------------------------------------------------------------------------
int main(int argc, char** argv) {

  // Generate some random seed points.
  cout << "Generating random points....";
  clock_t t0 = clock();
  const unsigned n = 100000;
  double low[3] = {0.0, 0.0, 0.0};
  double high[3] = {1.0, 1.0, 1.0};
  vector<double> points;
  for (unsigned i = 0; i != n; ++i) {
    points.push_back(low[0] + (high[0] - low[0])*random01());
    points.push_back(low[1] + (high[1] - low[1])*random01());
    points.push_back(low[2] + (high[2] - low[2])*random01());
  }
  CHECK(points.size() == 3*n);
  clock_t t1 = clock();
  cout << "required " << double(t1 - t0)/CLOCKS_PER_SEC << " seconds." << endl;

  // Get the hull.
  cout << "Generating convex hull... ";
  t0 = clock();
  const double tolerance = 1.0e-10;
  polytope::PLC<3, double> hull = polytope::convexHull_3d(points, low, tolerance);
  t1 = clock();
  cout << "required " << double(t1 - t0)/CLOCKS_PER_SEC << " seconds." << endl;

  // Check that the hull contains all the input vertices.
  cout << "Checking hull for containment... ";
  t0 = clock();
  unsigned badFacet = 0, nFacets = hull.facets.size();
  for (unsigned i = 0; i != n; ++i) {
    badFacet = convexContains(hull, &points.front(), &points[3*i], tolerance);
    if (badFacet < nFacets) cerr << "Failed containment for facet : " << badFacet << endl;
    CHECK(badFacet == nFacets);
  }
  t1 = clock();
  cout << "required " << double(t1 - t0)/CLOCKS_PER_SEC << " seconds." << endl;

  // for (unsigned k = 0; k != hull.facets.size(); ++k) {
  //   cerr << "Facet " << k << " : ";
  //   for (unsigned j = 0; j != hull.facets[k].size(); ++j) cerr << " " << hull.facets[k][j];
  //   cerr << " : ";
  //   for (unsigned j = 0; j != hull.facets[k].size(); ++j) cerr << " (" << points[3*hull.facets[k][j]] << " " << points[3*hull.facets[k][j] + 1] << " " << points[3*hull.facets[k][j] + 2] << ")";
  //   cerr << endl;
  // }

  cout << "PASS" << endl;
  return 0;
}
