#include "polytope_write_silo.h"
#include "polytope_c.h"
#include "polytope.hh"

#if HAVE_SILO

using namespace polytope;

extern "C"
{

void polytope_write_silo(polytope_tessellation_t* mesh, 
                         int num_node_fields,
                         char** node_field_names,
                         polytope_real_t** node_fields,
                         int num_edge_fields,
                         char** edge_edge_names,
                         polytope_real_t** edge_fields,
                         int num_face_fields,
                         char** face_field_names,
                         polytope_real_t** face_fields,
                         int num_cell_fields,
                         char** cell_field_names,
                         polytope_real_t** cell_fields,
                         const char* file_prefix,
                         const char* directory,
                         int cycle,
                         polytope_real_t time,
                         MPI_Comm comm,
                         int num_files,
                         int mpi_tag)
{
  // FIXME
}

}

#endif

