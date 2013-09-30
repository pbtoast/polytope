#include "polytope_tessellation.h"
#include "Tessellation.hh"

using namespace std;
using namespace polytope;

extern "C"
{

struct polytope_tessellation_t 
{
  Tessellation<2, polytope_real_t>* tess2;
  Tessellation<3, polytope_real_t>* tess3;
};

//------------------------------------------------------------------------
void polytope_tessellation_data_free(polytope_tessellation_data_t* tessellation_data)
{
  free(tessellation_data->nodes);
  free(tessellation_data->cell_offsets);
  free(tessellation_data->cell_faces);
  free(tessellation_data->face_offsets);
  free(tessellation_data->face_nodes);
  free(tessellation_data->inf_nodes);
  free(tessellation_data->inf_faces);
  free(tessellation_data->face_cells);
  if (tessellation_data->neighbor_domains != NULL)
    free(tessellation_data->neighbor_domains);
  if (tessellation_data->shared_node_domain_offsets != NULL)
  {
    free(tessellation_data->shared_node_domain_offsets);
    free(tessellation_data->shared_nodes);
  }
  if (tessellation_data->shared_face_domain_offsets != NULL)
  {
    free(tessellation_data->shared_face_domain_offsets);
    free(tessellation_data->shared_faces);
  }
  free(tessellation_data);
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_tessellation_fprintf(polytope_tessellation_t* tessellation, FILE* stream)
{
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_tessellation_free(polytope_tessellation_t* tessellation)
{
  if (tessellation->tess2 != NULL)
    delete tessellation->tess2;
  else
    delete tessellation->tess3;
  free(tessellation);
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
polytope_tessellation_data_t* polytope_tessellation_data(polytope_tessellation_t* tessellation)
{

}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
polytope_plc_t* polytope_tessellation_convex_hull(polytope_tessellation_t* tessellation)
{
  return NULL;
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_get_num_node_cells(polytope_tessellation_t* tessellation,
                                 int* num_node_cells)
{
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void polytope_get_node_cells(polytope_tessellation_t* tessellation,
                             int* node_cells)
{
}
//------------------------------------------------------------------------

}

