//------------------------------------------------------------------------------
// Handy method to make a box PLC from a low & high point.
//------------------------------------------------------------------------------
#ifndef __polytope_makeBoxPLC__
#define __polytope_makeBoxPLC__

#include "ReducedPLC.hh"

namespace polytope {
namespace internal {

//-------------------------------------------------------------------------------
// 2D
//-------------------------------------------------------------------------------
template<typename RealType>
void
makeBoxPLC(ReducedPLC<2, RealType>& result,
           const RealType* low, const RealType* high) {
  result.points.resize(8);
  result.points[0] =  low[0]; result.points[1] =  low[1];
  result.points[2] = high[0]; result.points[3] =  low[1];
  result.points[4] = high[0]; result.points[5] = high[1];
  result.points[6] =  low[0]; result.points[7] = high[1];
  result.facets.resize(4, std::vector<int>(2));
  for (unsigned i = 0; i != 4; ++i) {
    result.facets[i][0] = i;
    result.facets[i][1] = (i + 1) % 4;
  }
}

//-------------------------------------------------------------------------------
// 3D
//-------------------------------------------------------------------------------
template<typename RealType>
void
makeBoxPLC(ReducedPLC<3, RealType>& result,
           const RealType* low, const RealType* high) {
  POLY_ASSERT2(false, "Implement me!");
}

}
}

#endif
