#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "polytope_c.h"
#include "polytope_test_utilities.h"

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#define MIN(x, y) ((x < y) ? x : y)
#define MAX(x, y) ((x > y) ? x : y)

//------------------------------------------------------------------------
// test1
//------------------------------------------------------------------------
void test1(polytope_tessellator_t* tessellator) {
  
  // Output name
  const char* test_name = "boost_tessellator_test1";

  // Create the generators.
  const int nx = 3;
  int num_points = nx*nx;
  double points[2*num_points];
  const double x1 = 0.0, y1 = 0.0;
  const double x2 = 1.0, y2 = 1.0;
  const double dx = (x2 - x1)/nx, dy = (y2 - y1)/nx;
  double low [2] = { FLT_MAX, FLT_MAX};
  double high[2] = {-FLT_MAX, -FLT_MAX};
  unsigned ix, iy, offset = 0;
  double xi, yi;
  for (iy = 0; iy != nx; ++iy) {
    yi = y1 + (iy + 0.5)*dy;
    for (ix = 0; ix != nx; ++ix) {
      xi = x1 + (ix + 0.5)*dx;
      points[offset++] = xi;
      points[offset++] = yi;
      low [0] = MIN(low [0], xi);
      low [1] = MIN(low [1], yi);
      high[0] = MAX(high[0], xi);
      high[1] = MAX(high[1], yi);
    }
  }

  // Construct a PLC to represent the boundary.
  int num_plc_points = 8;
  double plc_points[8] = {x1, y1, x2, y1, x2, y2, x1, y2};
  polytope_plc_t* plc = polytope_plc_new(2);
  polytope_plc_add_facet(plc);
  polytope_plc_add_facet(plc);
  polytope_plc_add_facet(plc);
  polytope_plc_add_facet(plc);
  int i;
  for (i = 0; i < 4; ++i) 
  {
    polytope_plc_add_facet_node(plc, i, i);
    polytope_plc_add_facet_node(plc, i, (i+1)%4);
  }

  polytope_tessellation_t* mesh = polytope_tessellation_new(2);

  // Tessellate unbounded
  polytope_tessellator_tessellate_unbounded(tessellator, points, num_points, mesh);
  output_mesh(mesh, test_name, points, num_points, 1, 0.0);
  polytope_tessellation_clear(mesh);

  // Tessellate bounded
  polytope_tessellator_tessellate_in_plc(tessellator, points, num_points, 
                                         plc_points, num_plc_points,
                                         plc, mesh);
  output_mesh(mesh,test_name, points, num_points, 2, 0.0);

  polytope_plc_free(plc);
}


//------------------------------------------------------------------------
// test2
//------------------------------------------------------------------------
void test2(polytope_tessellator_t* tessellator) 
{
#if 0  
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
  cout << mesh << endl;
  mesh.clear();

  // Tessellate bounded
  tessellator.tessellate(points, PLCpoints, boundary, mesh);
  outputMesh(mesh,testName,points,2);

  //cout << mesh << endl;
#endif
}

//------------------------------------------------------------------------
// main
//------------------------------------------------------------------------
int 
main(int argc, char** argv) 
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  polytope_tessellator_t* tessellator = boost_tessellator_new();
  
  {
    printf("\nTest 1\n");
    test1(tessellator);
  }

  {
    printf("\nTest 2\n");
    test2(tessellator);
  }

  polytope_tessellator_free(tessellator);
  printf("PASS\n");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return 0;
}
//------------------------------------------------------------------------
