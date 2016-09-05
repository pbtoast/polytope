//------------------------------------------------------------------------------
// Clip a QuantizedTessellation against a PLC boundary.
//------------------------------------------------------------------------------
#ifndef __polytope_clipQuantizedTessellation__
#define __polytope_clipQuantizedTessellation__

#include "QuantizedTessellation2d.hh"
#include "QuantizedTessellation3d.hh"
#include "PLC.hh"

namespace polytope {

// 2D
template<typename IntType, typename RealType>
void clipQuantizedTessellation(QuantizedTessellation2d<IntType, RealType>& qmesh,
                               const std::vector<RealType>& PLCpoints,
                               const PLC<2, RealType>& geometry);

// 3D
template<typename IntType, typename RealType>
void clipQuantizedTessellation(QuantizedTessellation3d<IntType, RealType>& qmesh,
                               const std::vector<RealType>& PLCpoints,
                               const PLC<3, RealType>& geometry);

}

#endif
