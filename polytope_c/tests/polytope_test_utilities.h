//------------------------------------------------------------------------------
// C versions of random stuff useful for testing in polytope.
//------------------------------------------------------------------------------
#ifndef __polytope_c_test_utilities__
#define __polytope_c_test_utilities__

#include <float.h>
#include "polytope_c.h"

//------------------------------------------------------------------------------
// A macro for checking true/false test conditions.
//------------------------------------------------------------------------------
#define POLY_CHECK(x) if (!(x)) { sprintf("FAIL: %s\n", #x); exit(-1); }
#define POLY_CHECK2(x, msg) if (!(x)) { sprintf("FAIL: %s\n%s\n", #x, msg); exit(-1); }

//------------------------------------------------------------------------------
// Define our own local random number generator wrapping the standard srand &
// rand methods.
//------------------------------------------------------------------------------
double random01() {
  return 1.0*rand()/RAND_MAX;
}

//------------------------------------------------------------------------------
// A simple mesh output function for the SiloWriter
//------------------------------------------------------------------------------
void output_mesh(polytope_tessellation_t* mesh,
                 const char* prefix,
                 polytope_real_t* points,
                 int num_points,
                 unsigned test_cycle,
                 polytope_real_t time) 
{
#ifdef HAVE_SILO
  double* index = (double*)malloc(sizeof(double) * mesh->num_cells);
  double* genx = (double*)malloc(sizeof(double) * mesh->num_cells);
  double* geny = (double*)malloc(sizeof(double) * mesh->num_cells);
  int i;
  for (i = 0; i < mesh->num_cells; ++i)
  {
    index[i] = 1.0 * i;
    if (mesh->num_nodes != 0)
    {
      genx[i] = points[2*i  ];
      geny[i] = points[2*i+1];
    }
  }
  char *cell_field_names[mesh->dimension + 1];
  double *cell_fields[mesh->dimension + 1];
  cell_field_names[0] = (char*)"cell_index";
  cell_fields[0] = index;
  cell_field_names[1] = (char*)"gen_x";
  cell_fields[1] = genx;
  cell_field_names[2] = (char*)"gen_x";
  cell_fields[2] = genx;
  if (mesh->dimension == 3)
  {
    cell_field_names[3] = (char*)"gen_x";
    cell_fields[3] = genx;
  }
  polytope_write_silo(mesh, 0, NULL, NULL, 0, NULL, NULL, 0, NULL, NULL, mesh->dimension+1, cell_field_names, cell_fields,
                      prefix, ".", test_cycle, time, MPI_COMM_WORLD, -1, 0);
#endif
}

#endif
