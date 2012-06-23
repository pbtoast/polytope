// Try tessellating a simple lattice of generators in a box.

#include <iostream>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>

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
  triangle.tessellate(generators, generators, box, mesh);
  CHECK(mesh.nodes.size()/2 == (nx + 1)*(nx + 1));
  CHECK(mesh.cells.size() == nx*nx);
  for (unsigned i = 0; i != nx*nx; ++i) 
  {
    CHECK(mesh.cells[i].size() == 4);
  }
  CHECK(mesh.faces.size() == 2*nx*(nx + 1));

  // Write out the file if we can.
#ifdef HAVE_SILO
  vector<double> index(mesh.cells.size());
  for (int i = 0; i < mesh.cells.size(); ++i) index[i] = double(i);
  map<string, double*> fields;
  fields["cell_index"] = &index[0];
  SiloWriter<2, double>::write(mesh, fields, "test_TriangleTessellator_2x2");
#endif
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void
testCircle()
{
  // Generate a bunch of random points within a circle of radius 1.
  vector<double> generators, PLCpoints;
  int N = 100;
  for (int i = 0; i < N; ++i)
  {
    double x = FLT_MAX, y = FLT_MAX;
    double safetyBuffer = 0.05; // Important!!
    while (x*x + y*y >= 1.0 - safetyBuffer)
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
    PLCpoints.push_back(x);
    PLCpoints.push_back(y);
  }

  // Facets.
  PLC<2, double> circle;
  circle.facets.resize(Nb); 
  for (int f = 0; f < Nb - 1; ++f)
  {
    circle.facets[f].resize(2);
    circle.facets[f][0] = f;
    circle.facets[f][1] = f+1;
  }

  // Last facet completes the circle.
  int f = Nb-1;
  circle.facets[f].resize(2);
  circle.facets[f][0] = f;
  circle.facets[f][1] = 0;

  // Create the tessellation.
  Tessellation<2, double> mesh;
  TriangleTessellator<double> triangle;
  triangle.tessellate(generators, PLCpoints, circle, mesh);
  CHECK(mesh.cells.size() == N);

  // // Make sure that the nodes all fall within the circle.
  // for (int n = 0; n < mesh.nodes.size()/2; ++n)
  // {
  //   double x = mesh.nodes[2*n], y = mesh.nodes[2*n+1];
  //   CHECK(x*x + y*y < 1+1e-14);
  // }

  // Write out the file if we can.
#ifdef HAVE_SILO
  vector<double> index(mesh.cells.size());
  for (int i = 0; i < mesh.cells.size(); ++i) index[i] = double(i);
  map<string, double*> fields;
  fields["cell_index"] = &index[0];
  SiloWriter<2, double>::write(mesh, fields, "test_TriangleTessellator_circle");
#endif
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void
testDonut()
{
  // Generate a bunch of random points within a circle of radius 1.
  vector<double> generators, PLCpoints;
  int N = 100;
  for (int i = 0; i < N; ++i)
  {
    double x = FLT_MAX, y = FLT_MAX;
    double outerBuffer = 0.05, innerBuffer = 0.01; // Important!!
    while ((x*x + y*y >= (1.0 - outerBuffer)*(1.0 - outerBuffer)) or
           (x*x + y*y <= (0.25 + innerBuffer)*(0.25 + innerBuffer)))
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

  // Outer circle.
  for (int b = 0; b < Nb; ++b)
  {
    double theta = 2.0*M_PI*double(b)/double(Nb+1);
    double x = cos(theta), y = sin(theta);
    PLCpoints.push_back(x);
    PLCpoints.push_back(y);
  }

  // Facets on the outer circle.
  PLC<2, double> donut;
  donut.facets.resize(Nb); 
  for (int f = 0; f < Nb - 1; ++f)
  {
    donut.facets[f].resize(2);
    donut.facets[f][0] = PLCpoints.size()/2-Nb+f;
    donut.facets[f][1] = PLCpoints.size()/2-Nb+f+1;
  }

  // Last facet completes the outer circle.
  int f = Nb-1;
  donut.facets[f].resize(2);
  donut.facets[f][0] = PLCpoints.size()/2-Nb+f;
  donut.facets[f][1] = PLCpoints.size()/2-Nb;

  // Inner circle.
  for (int b = 0; b < Nb; ++b)
  {
    double theta = 2.0*M_PI*(1.0 - double(b)/double(Nb+1));
    double x = 0.25*cos(theta), y = 0.25*sin(theta);
    PLCpoints.push_back(x);
    PLCpoints.push_back(y);
  }

  // Facets on the inner circle.
  donut.holes = vector<vector<vector<int> > >(1);
  donut.holes[0].resize(Nb);
  for (int f = 0; f < Nb - 1; ++f)
  {
    donut.holes[0][f].resize(2);
    donut.holes[0][f][0] = PLCpoints.size()/2-Nb+f;
    donut.holes[0][f][1] = PLCpoints.size()/2-Nb+f+1;
  }

  // Last facet completes the inner circle.
  f = Nb-1;
  donut.holes[0][f].resize(2);
  donut.holes[0][f][0] = PLCpoints.size()/2-Nb+f;
  donut.holes[0][f][1] = PLCpoints.size()/2-Nb;

  // Create the tessellation.
  Tessellation<2, double> mesh;
  TriangleTessellator<double> triangle;
  triangle.tessellate(generators, PLCpoints, donut, mesh);

//   // Make sure that the nodes all fall within the donut.
//   for (int n = 0; n < mesh.nodes.size()/2; ++n)
//   {
//     double x = mesh.nodes[2*n], y = mesh.nodes[2*n+1];
//     CHECK(x*x + y*y < 1+1e-14);
// //    CHECK(x*x + y*y > 0.25*0.25-1e-14);
//   }

  // Write out the file if we can.
#ifdef HAVE_SILO
  vector<double> index(mesh.cells.size());
  for (int i = 0; i < mesh.cells.size(); ++i) index[i] = double(i);
  map<string, double*> fields;
  fields["cell_index"] = &index[0];
  SiloWriter<2, double>::write(mesh, fields, "test_TriangleTessellator_donut");
#endif
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void
testBounded()
{
  // Generate a bunch of random points.
  vector<double> generators;
  int N = 100;
  double L = 10.0;
  set<int> used;
  const int ixmax = 1 << 14;
  const int imax = ixmax*ixmax;
  ASSERT(imax < RAND_MAX);
  int ix, iy, iran;
  double x, y;
  for (int i = 0; i < N; ++i)
  {
    iran = ::random() % imax;
    while (used.find(iran) != used.end()) iran = ::random() % imax;
    used.insert(iran);
    ix = iran % ixmax;
    iy = iran / ixmax;
    x = L * (double(ix)/ixmax - 0.5);
    y = L * (double(iy)/ixmax - 0.5);
    generators.push_back(x);
    generators.push_back(y);
  }

  // Create a bounding box for the square.
  double low[2] = {-0.55*L, -0.55*L},
        high[2] = { 0.55*L,  0.55*L};

  // Create the tessellation.
  Tessellation<2, double> mesh;
  TriangleTessellator<double> triangle;
  triangle.tessellate(generators, low, high, mesh);
  CHECK(mesh.cells.size() == generators.size()/2);

  // Write out the file if we can.
#ifdef HAVE_SILO
  vector<double> index(mesh.cells.size());
  for (int i = 0; i < mesh.cells.size(); ++i) index[i] = double(i);
  map<string, double*> fields;
  fields["cell_index"] = &index[0];
  SiloWriter<2, double>::write(mesh, fields, "test_TriangleTessellator_bounded");
#endif
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void
testUnbounded()
{
  // Generate a bunch of random points.
  vector<double> generators;
  int N = 4;
  double L = 10.0;
  set<int> used;
  const int ixmax = 1 << 14;
  const int imax = ixmax*ixmax;
  ASSERT(imax < RAND_MAX);
  int ix, iy, iran;
  double x, y;
  for (int i = 0; i < N; ++i)
  {
    iran = ::random() % imax;
    while (used.find(iran) != used.end()) iran = ::random() % imax;
    used.insert(iran);
    ix = iran % ixmax;
    iy = iran / ixmax;
    x = L * (double(ix)/ixmax - 0.5);
    y = L * (double(iy)/ixmax - 0.5);
    generators.push_back(x);
    generators.push_back(y);
  }

  // Blago!
  generators.resize(0);
  generators.push_back(-1.0); generators.push_back(0.0);
  generators.push_back( 1.0); generators.push_back(0.0);
  generators.push_back( 0.0); generators.push_back(0.2);
  generators.push_back( 0.9); generators.push_back(0.1);
  // Blago!

  // Create the tessellation.
  Tessellation<2, double> mesh;
  TriangleTessellator<double> triangle;
  triangle.tessellate(generators, mesh);
  CHECK(mesh.cells.size() == N);

  // Write out the file if we can.
#ifdef HAVE_SILO
  vector<double> index(mesh.cells.size());
  for (int i = 0; i < mesh.cells.size(); ++i) index[i] = double(i);
  map<string, double*> fields;
  fields["cell_index"] = &index[0];
  SiloWriter<2, double>::write(mesh, fields, "test_TriangleTessellator_unbounded");
#endif
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
int 
main(int argc, char** argv) 
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  test2x2Box();
  testCircle();
  testDonut();
  testBounded(); 
  testUnbounded();

  cout << "PASS" << endl;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
//------------------------------------------------------------------------
