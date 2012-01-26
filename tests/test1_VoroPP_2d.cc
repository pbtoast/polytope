// Try tessellating a simple lattice of generators in a box.

#include <iostream>
#include <vector>
#include <cassert>

#include "VoroPP_2d.hh"

using namespace std;

int main() {

  // Create the generators.
  vector<double> generators;
  const double nx = 2;
  const double x1 = 0.0, y1 = 0.0;
  const double x2 = 100.0, y2 = 100.0;
  const double dx = (x2 - x1)/nx, dy = (y2 - y1)/nx;
  unsigned ix, iy;
  double xi, yi;
  for (iy = 0; iy != nx; ++iy) {
    yi = y1 + (iy + 0.5)*dx;
    for (ix = 0; ix != nx; ++ix) {
      xi = x1 + (ix + 0.5)*dx;
      generators.push_back(xi);
      generators.push_back(yi);
    }
  }

  // Create the tessellation.
  polytope::Tessellation<double> mesh;
  polytope::VoroPP_2d<double> voro(x1, y1, x2, y2);
  voro.tessellate(generators, mesh);
  assert(mesh.nodes.size() == (nx + 1)*(nx + 1));
  assert(mesh.cells.size() == nx*nx);
  for (unsigned i = 0; i != nx*nx; ++i) assert(mesh.cells[i].size() == 4);
  assert(mesh.faces.size() == 2*nx*(nx + 1));

  cout << "PASS" << endl;
  return 0;
}
