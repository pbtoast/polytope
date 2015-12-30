#ifndef POLYTOPE_C_READ_SILO_H
#define POLYTOPE_C_READ_SILO_H

#ifdef HAVE_SILO

#ifdef HAVE_MPI
#include <mpi.h>
#else
#ifndef MPI_Comm
#define MPI_Comm int
#define MPI_COMM_WORLD 0
#endif
#endif 

#ifdef __cplusplus
extern "C"
{
#endif

//! Reads an arbitrary polygonal mesh and an associated set of 
//! cell-centered fields to a SILO file in the given directory.
void polytope_read_silo(polytope_tessellation_t* mesh, 
                        int* num_fields,
                        char*** field_names,
                        polytope_real_t*** fields,
                        const char* file_prefix,
                        const char* directory,
                        int cycle,
                        polytope_real_t* time,
                        MPI_Comm comm,
                        int num_files,
                        int mpi_tag);

//! This version allows tags to be read in addition to the mesh and field data.
void polytope_read_silo_with_tags(polytope_tessellation_t* mesh, 
                                  int* num_fields,
                                  char*** field_names,
                                  polytope_real_t*** fields,
                                  int* num_node_tags,
                                  char*** node_tag_names,
                                  int** node_tag_sizes,
                                  int*** node_tags,
                                  int* num_edge_tags,
                                  char*** edge_tag_names,
                                  int** edge_tag_sizes,
                                  int*** edge_tags,
                                  int* num_face_tags,
                                  char*** face_tag_names,
                                  int** face_tag_sizes,
                                  int*** face_tags,
                                  int* num_cell_tags,
                                  char*** cell_tag_names,
                                  int** cell_tag_sizes,
                                  int*** cell_tags,
                                  const char* file_prefix,
                                  const char* directory,
                                  int cycle,
                                  polytope_real_t* time,
                                  MPI_Comm comm,
                                  int num_files,
                                  int mpi_tag);

#ifdef __cplusplus
}
#endif

#endif

#endif
