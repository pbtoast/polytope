//-----------------------------------------------------------------------------//
// PLC_CSG instantiation and declarations.
//----------------------------------------------------------------------------//
#include "polytope.hh"
#include "PLC_CSG_2d.hh"
#include "PLC_CSG_3d.hh"
#include "polytope_geometric_utilities.hh"

namespace polytope {
  namespace CSG {
    namespace CSG_internal_2d {
      template<> double Line<float>::EPSILON = 1.0/geometry::Hasher<2, float>::coordMax();
      template<> double Line<double>::EPSILON = 1.0/geometry::Hasher<2, double>::coordMax();
      template<> double Line<int64_t>::EPSILON = 1.0/geometry::Hasher<2, int64_t>::coordMax();
    }
    namespace CSG_internal_3d {
      template<> double Plane<float>::EPSILON = 1.0/geometry::Hasher<3, float>::coordMax();
      template<> double Plane<double>::EPSILON = 1.0/geometry::Hasher<3, double>::coordMax();
      template<> double Plane<int64_t>::EPSILON = 1.0/geometry::Hasher<3, int64_t>::coordMax();
    }
  }
}
