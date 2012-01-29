// Try tessellating a simple lattice of generators in a box.

#include <iostream>
#include <vector>
#include <cassert>

#include "polytope.hh"

#define CHECK(x) if (!(x)) { cout << "FAIL: " << #x << endl; exit(-1); }

using namespace std;
using namespace polytope;

//------------------------------------------------------------------------
void
test2x2Box()
{
  // Create the generators.
  vector<double> generators;
  const int nx = 2;
  const double x1 = 0.0, y1 = 0.0;
  const double x2 = 100.0, y2 = 100.0;
  const double dx = (x2 - x1)/nx, dy = (y2 - y1)/nx;
  unsigned ix, iy;
  double xi, yi;
  for (iy = 0; iy != nx; ++iy) {
    yi = y1 + (iy + 0.5)*dy;
    for (ix = 0; ix != nx; ++ix) {
      xi = x1 + (ix + 0.5)*dx;
      generators.push_back(xi);
      generators.push_back(yi);
    }
  }

  // Create the piecewise linear complex representing the box. Note that 
  // the box consists of facets that are defined by their connections to 
  // generating points.
  PLC<2, double> box;

  // 4 facets
  box.facets.resize(4);

  // facet 0 -- connects 0th and (nx-1) generator points.
  box.facets[0].resize(2);
  box.facets[0][0] = 0; box.facets[0][1] = nx-1;

  // facet 1 -- connects (nx-1)th and (nx*nx-1)th generator points.
  box.facets[1].resize(2);
  box.facets[1][0] = nx-1; box.facets[1][1] = nx*nx-1;

  // facet 2 -- connects (nx*nx-1)th and (nx*(nx-1))th generator points.
  box.facets[2].resize(2);
  box.facets[2][0] = nx*nx-1; box.facets[2][1] = nx*(nx-1);

  // facet 3 -- connects (nx*(nx-1))th and 0th generator points.
  box.facets[3].resize(2);
  box.facets[3][0] = nx*(nx-1); box.facets[3][1] = 0;

  // Create the tessellation.
  Tessellation<2, double> mesh;
  TriangleTessellator<double> triangle;
  triangle.tessellate(generators, box, mesh);
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
  SiloWriter<2, double>::write(mesh, fields, "test_TriangleTessellator_2x2");
#endif
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void
testCircle()
{
  // Generate a bunch of random points within a circle of radius 1.
  vector<double> generators;
  int N = 100;
  for (int i = 0; i < N; ++i)
  {
    double x = FLT_MAX, y = FLT_MAX;
    while (x*x + y*y >= 1.0)
    {
      x = 2.0 * double(::random())/RAND_MAX - 1.0;
      y = 2.0 * double(::random())/RAND_MAX - 1.0;
    }
    generators.push_back(x);
    generators.push_back(y);
  }

  // Create a piecewise linear complex representing the circle. Note that 
  // the PLC consists of facets that are defined by their connections to 
  // generating points, which means we need generators along the circle 
  // itself.

  // Boundary generators.
  int Nb = 90; // 4-degree resolution.
  for (int b = 0; b < Nb; ++b)
  {
    double theta = 2.0*M_PI*b/(Nb+1);
    double x = cos(theta), y = sin(theta);
    generators.push_back(x);
    generators.push_back(y);
  }

  // Facets.
  PLC<2, double> circle;
  circle.facets.resize(Nb); 
  for (int f = 0; f < Nb - 1; ++f)
  {
    circle.facets[f].resize(2);
    circle.facets[f][0] = generators.size()/2-Nb+f;
    circle.facets[f][1] = generators.size()/2-Nb+f+1;
  }

  // Last facet completes the circle.
  int f = Nb-1;
  circle.facets[f].resize(2);
  circle.facets[f][0] = generators.size()/2-Nb+f;
  circle.facets[f][1] = generators.size()/2-Nb;

  // Create the tessellation.
  Tessellation<2, double> mesh;
  TriangleTessellator<double> triangle;
  triangle.tessellate(generators, circle, mesh);
  CHECK(mesh.cells.size() == (N + Nb));

  // Make sure that the nodes all fall within the circle.
  for (int n = 0; n < mesh.nodes.size()/2; ++n)
  {
    double x = mesh.nodes[2*n], y = mesh.nodes[2*n+1];
    CHECK(x*x + y*y < 1);
  }

//  for (unsigned i = 0; i != nx*nx; ++i) 
//  {
//    CHECK(mesh.cells[i].size() == 4);
//  }
//  CHECK(mesh.faces.size() == 2*nx*(nx + 1));

  // Write out the file if we can.
#ifdef HAVE_SILO
  vector<double> r2(N+Nb);
  for (int i = 0; i < N+Nb; ++i)
    r2[i] = double(::random())/RAND_MAX;
  map<string, double*> fields;
  fields["data"] = &r2[0];
  SiloWriter<2, double>::write(mesh, fields, "test_TriangleTessellator_circle");
#endif
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
int 
main() 
{
  test2x2Box();
  testCircle();

  cout << "PASS" << endl;
  return 0;
}
//------------------------------------------------------------------------
