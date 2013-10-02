#ifndef POLYTOPE_C_WRITE_SILO_H
#define POLYTOPE_C_WRITE_SILO_H

#if HAVE_SILO

#if HAVE_MPI
#include <mpi.h>
#else
#ifndef MPI_Comm
#define MPI_Comm int
#endif
#endif 

#ifdef __cplusplus
extern "C"
{
#endif

// Writes an arbitrary polyhedral mesh and an associated set of 
// (node, edge, face, cell)-centered fields to a SILO file in the given directory.
void polytope_write_silo(polytope_tessellation_t* mesh, 
                         int num_node_fields,
                         char** node_field_names,
                         polytope_real_t** node_fields,
                         int num_edge_fields,
                         char** edge_field_names,
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
                         int mpi_tag);

#ifdef __cplusplus
}
#endif

#endif

#endif
