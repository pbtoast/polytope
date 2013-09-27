#include "polytope_read_silo.h"
#include "polytope_c.h"
#include "polytope.hh"

#if HAVE_SILO

extern "C"
{

void polytope_read_silo(polytope_tessellation_t* mesh, 
                        int* num_fields,
                        char*** field_names,
                        polytope_real_t*** fields,
                        const char* file_prefix,
                        const char* directory,
                        int cycle,
                        RealType& time,
                        MPI_Comm comm,
                        int num_files,
                        int mpi_tag)
{
  // FIXME
}

}

#endif

