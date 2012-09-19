#include <mpi.h>

namespace polytope {

void internal_abort() {
  MPI_Abort(MPI_COMM_WORLD, -1);

}
}

