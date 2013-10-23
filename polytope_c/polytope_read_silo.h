#ifndef POLYTOPE_C_READ_SILO_H
#define POLYTOPE_C_READ_SILO_H

#if HAVE_SILO

#if HAVE_MPI
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
                        polytope_real_t time,
                        MPI_Comm comm,
                        int num_files,
                        int mpi_tag);

#ifdef __cplusplus
}
#endif

#endif

#endif
