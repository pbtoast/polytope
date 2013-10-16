//-----------------------------------------------------------------------------//
// PLC_CSG instantiation and declarations.
//----------------------------------------------------------------------------//
#include "polytope.hh"
#include "PLC_CSG_3d.hh"

namespace polytope {
  namespace CSG {
    namespace CSG_internal_3d {
      template<> double Plane<float>::EPSILON = 1.0e-8;
      template<> double Plane<double>::EPSILON = 1.0e-10;
      template<> double Plane<int64_t>::EPSILON = 1.0e-8;
    }
  }
}
