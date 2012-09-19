#include <iostream>
#include <cstdlib>

#include "polytope_internal.hh"

#if HAVE_MPI
#include <mpi.h>
#endif

namespace polytope {

void internal_abort() {
  std::cout.flush();
  std::cerr.flush();

#if HAVE_MPI
  MPI_Abort(MPI_COMM_WORLD, -1);
#else
  abort();
#endif
}

}
