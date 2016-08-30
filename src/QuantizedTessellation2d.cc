//------------------------------------------------------------------------
// QuantizedTessellation2d
//
// An intermediate representation for 2D tessellations in integer
// coordinates.
//------------------------------------------------------------------------
#include <iostream>
#include <algorithm>
#include <map>

#include "QuantizedTessellation2d.hh"
#include "polytope.hh"
#include "polytope_tessellator_utilities.hh"

namespace polytope {
  
//------------------------------------------------------------------------------
// Construct with the given generators.  Finds the bounding limits, sets the 
// infRadius, and sets the quantized generators.
//------------------------------------------------------------------------------
template<typename IntType, typename RealType>
QuantizedTessellation2d<IntType, RealType>::
QuantizedTessellation2d(const std::vector<RealType>& points) {
  geometry::computeBoundingBox<2, RealType>(&points[0], points.size(), true, xmin, xmax);
  length = std::max(xmax[0] - xmin[0], xmax[1] - xmin[1]);
  xmin[0] -= 2.0*length;
  xmin[1] -= 2.0*length;
  xmax[0] += 2.0*length;
  xmax[1] += 2.0*length;
  infRadius = 1.5*length;
  length *= 5.0;
  const int numGenerators = points.size()/2;
  generators.resize(numGenerators);
  for (unsigned i = 0; i < numGenerators; ++i) {
    this->quantize(&points[2*i], &generators[i].x);
    generators[i].index = i;
  }
}

//------------------------------------------------------------------------------
// Convert real coordinates to integers.
//------------------------------------------------------------------------------
template<typename IntType, typename RealType>
void
QuantizedTessellation2d<IntType, RealType>::
quantize(const RealType* realcoords, IntType* intcoords) const {
  const RealType dx = length/(std::numeric_limits<IntType>::max()/4 - std::numeric_limits<IntType>::min()/4);
  intcoords[0] = std::numeric_limits<IntType>::min()/4 + IntType((realcoords[0] - xmin[0])/dx);
  intcoords[1] = std::numeric_limits<IntType>::min()/4 + IntType((realcoords[1] - xmin[1])/dx);
}

//------------------------------------------------------------------------------
// Convert int coordinates to reals.
//------------------------------------------------------------------------------
template<typename IntType, typename RealType>
void
QuantizedTessellation2d<IntType, RealType>::
dequantize(const IntType* intcoords, RealType* realcoords) {
  const RealType dx = length/(std::numeric_limits<IntType>::max()/4 - std::numeric_limits<IntType>::min()/4);
  realcoords[0] = xmin[0] + (intcoords[0] - std::numeric_limits<IntType>::min()/4)*dx;
  realcoords[1] = xmin[1] + (intcoords[1] - std::numeric_limits<IntType>::min()/4)*dx;
}

//------------------------------------------------------------------------------
// Instantiate versions we know we need.
//------------------------------------------------------------------------------
template struct QuantizedTessellation2d<int, double>;

}
