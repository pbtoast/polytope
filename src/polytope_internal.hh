// polytope_internal.hh
//
// Put common includes for polytope here that you don't necessarily 
// want exposed in the public interface.

// An ASSERT macro, if one isn't already defined.
#ifndef ASSERT

#if HAVE_MPI

// Parallel ASSERT -- calls MPI_Abort if NDEBUG is defined.
#include <mpi.h>
#include <iostream>

// Forward declare our helper abort method.
namespace polytope {
  void internal_abort();
}

#ifndef NDEBUG
#define ASSERT(x) \
  if (!(x)) \
  { \
    std::cout << "Assertion " << #x << " failed\nat " << __FILE__ << ":" << __LINE__ << std::endl; \
    polytope::internal_abort(); \
  }
#else
#define ASSERT(x) 
#endif

#else

// Serial ASSERT -- forwarded to C assert.
#include <cassert>
#define ASSERT(x) assert(x)
#endif

#endif

// Classes within the library.
#include "ErrorHandler.hh"
