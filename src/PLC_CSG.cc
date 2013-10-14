//-----------------------------------------------------------------------------//
// PLC_CSG instantiation and declarations.
//----------------------------------------------------------------------------//
#include "polytope.hh"
#include "PLC_CSG.hh"

namespace polytope {
  namespace CSG {
    namespace CSG_internal {
      template<> float  Plane<float>::EPSILON  = 1.0e-8f;
      template<> double Plane<double>::EPSILON = 1.0e-10;
      template<> int64_t Plane<int64_t>::EPSILON = 2;
    }
  }
}
