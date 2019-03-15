#include <stdlib.h>
#include <string.h>

#include "polytope_c.h"
#include "Tessellation.hh"

using namespace std;
using namespace polytope;

extern "C"
{

//------------------------------------------------------------------------
polytope_tessellation_t* polytope_tessellation_new(int dimension)
{
  POLY_ASSERT((dimension == 2) or (dimension == 3));
  polytope_tessellation_t* tess = (polytope_tessellation_t*)malloc(sizeof(polytope_tessellation_t));
  tess->dimension = dimension;
  tess->num_nodes = 0;
  tess->nodes = NULL;
  tess->num_cells = 0;
  tess->cell_offsets = NULL;
  tess->cell_faces = NULL;
  tess->num_faces = 0;
  tess->face_offsets = NULL;
  tess->face_nodes = NULL;
  tess->num_neighbor_domains = 0;
  tess->neighbor_domains = NULL;
  tess->shared_node_domain_offsets = NULL;
  tess->shared_nodes = NULL;
  tess->shared_face_domain_offsets = NULL;
  tess->shared_faces = NULL;
  tess->node_cell_offsets = NULL;
  tess->node_cells = NULL;
  tess->convex_hull = NULL;
  tess->boundary_faces = NULL;
  tess->num_boundary_faces = 0;
  tess->boundary_nodes = NULL;
  tess->num_boundary_nodes = 0;

  return tess;
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_tessellation_free(polytope_tessellation_t* tessellation)
{
  polytope_tessellation_clear(tessellation);
  free(tessellation);
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_tessellation_clear(polytope_tessellation_t* tessellation)
{
  if (tessellation->nodes != NULL)
  {
    free(tessellation->nodes);
    tessellation->num_nodes = 0;
    tessellation->nodes = NULL;
  }
  if (tessellation->cell_faces != NULL)
  {
    free(tessellation->cell_offsets);
    free(tessellation->cell_faces);
    tessellation->num_cells = 0;
    tessellation->cell_offsets = NULL;
    tessellation->cell_faces = NULL;
  }
  if (tessellation->face_nodes != NULL)
  {
    free(tessellation->face_offsets);
    free(tessellation->face_nodes);
    tessellation->num_faces = 0;
    tessellation->face_offsets = NULL;
    tessellation->face_nodes = NULL;
  }
  if (tessellation->boundary_nodes != NULL)
  {
    free(tessellation->boundary_nodes);
  }
  if (tessellation->boundary_faces != NULL)
  {
    free(tessellation->boundary_faces);
  }
  if (tessellation->face_cells != NULL)
  {
    free(tessellation->face_cells);
  }
  if (tessellation->neighbor_domains != NULL)
  {
    free(tessellation->neighbor_domains);
    tessellation->num_neighbor_domains = 0;
    tessellation->neighbor_domains = NULL;
  }
  if (tessellation->shared_node_domain_offsets != NULL)
  {
    free(tessellation->shared_node_domain_offsets);
    free(tessellation->shared_nodes);
    tessellation->shared_node_domain_offsets = NULL;
    tessellation->shared_nodes = NULL;
  }
  if (tessellation->shared_face_domain_offsets != NULL)
  {
    free(tessellation->shared_face_domain_offsets);
    free(tessellation->shared_faces);
    tessellation->shared_face_domain_offsets = NULL;
    tessellation->shared_faces = NULL;
  }
  if (tessellation->node_cell_offsets != NULL)
  {
    free(tessellation->node_cell_offsets);
    free(tessellation->node_cells);
  }
  if (tessellation->convex_hull != NULL)
  {
    polytope_plc_free(tessellation->convex_hull);
    tessellation->convex_hull = NULL;
  }

  // FIXME: Not used yet.
  tessellation->node_cell_offsets = NULL;
  tessellation->node_cells = NULL;
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_tessellation_snprintf(polytope_tessellation_t* tessellation, char* str, int n)
{
  POLY_ASSERT(tessellation != NULL);

  char fmt[128];
  snprintf(fmt, 128, "Tessellation (%d)", static_cast<int>(sizeof(str) - strlen(str) - 1));
  strncat(str, fmt, sizeof(str) - strlen(str));
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_tessellation_fprintf(polytope_tessellation_t* tessellation, FILE* stream)
{
  POLY_ASSERT(tessellation != NULL);
  // Let's assume, for the moment, that the output of our tessellation scales with the
  // number of cells.
  int n = 512 * tessellation->num_cells;
  char* str = (char*)malloc(n * sizeof(char));
  polytope_tessellation_snprintf(tessellation, str, n);
  fprintf(stream, "%s", str);
  free(str);
}
//------------------------------------------------------------------------

}
