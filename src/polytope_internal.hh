// polytope_internal.hh
//
// Put common includes for polytope here that you don't necessarily 
// want exposed in the public interface.

// An ASSERT macro, if one isn't already defined.
#ifndef ASSERT

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
#define ASSERT2(x, msg) \
  if (!(x)) \
  { \
    std::cout << "Assertion " << #x << " failed\nat " << __FILE__ << ":" << __LINE__ << std::endl << msg << std::endl; \
    polytope::internal_abort(); \
  }
#else
#define ASSERT(x)
#define ASSERT2(x, msg)
#endif

// Classes within the library.
#include "ErrorHandler.hh"

#endif
