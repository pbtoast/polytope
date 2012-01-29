// Try tessellating a simple lattice of generators in a box.

#include <iostream>
#include <vector>
#include <cassert>

#include "polytope.hh"

#define CHECK(x) if (!(x)) { cout << "FAIL: " << #x << endl; exit(-1); }

using namespace std;
using namespace polytope;

int main() {

  // Create the generators.
  vector<double> generators;
  const int nx = 2;
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
  Tessellation<double> mesh;
  TriangleTessellator<double> triangle;
  triangle.tessellate(generators, mesh);
  CHECK(mesh.nodes.size()/2 == (nx + 1)*(nx + 1));
  CHECK(mesh.cells.size() == nx*nx);
  for (unsigned i = 0; i != nx*nx; ++i) 
  {
    CHECK(mesh.cells[i].size() == 4);
  }
  CHECK(mesh.faces.size() == 2*nx*(nx + 1));

  // Write out the file if we can.
#ifdef HAVE_SILO
  vector<double> r2(nx*nx, 1.0);
  map<string, double*> fields;
  fields["data"] = &r2[0];
  SiloWriter<2, double>::write(mesh, fields, "test_TriangleTessellator");
#endif

  cout << "PASS" << endl;
  return 0;
}
