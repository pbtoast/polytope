//------------------------------------------------------------------------
// QuantizedTessellation3d
//
// An intermediate representation for 3D tessellations in integer
// coordinates.
//------------------------------------------------------------------------
#include <iostream>
#include <algorithm>
#include <map>

#include "QuantizedTessellation3d.hh"
#include "polytope.hh"
#include "polytope_tessellator_utilities.hh"

namespace polytope {
  
//------------------------------------------------------------------------------
// Construct with the given generators.  Finds the bounding limits, sets the 
// infRadius, and sets the quantized generators.
//------------------------------------------------------------------------------
template<typename IntType, typename RealType>
QuantizedTessellation3d<IntType, RealType>::
QuantizedTessellation3d(const std::vector<RealType>& points,
                        const std::vector<RealType>& boundaryPoints) {
  geometry::computeBoundingBox<3, RealType>(&boundaryPoints[0], boundaryPoints.size(), true, xmin, xmax);
  length = std::max(xmax[0] - xmin[0], std::max(xmax[1] - xmin[1], xmax[2] - xmin[2]));
  xmin[0] -= 2.0*length;
  xmin[1] -= 2.0*length;
  xmin[2] -= 2.0*length;
  xmax[0] += 2.0*length;
  xmax[1] += 2.0*length;
  xmax[2] += 2.0*length;
  infRadius = 1.5*length;
  length *= 5.0;
  const int numGenerators = points.size()/3;
  generators.resize(numGenerators);
  for (unsigned i = 0; i < numGenerators; ++i) {
    this->quantize(&points[3*i], &generators[i].x);
    generators[i].index = i;
  }
}

//------------------------------------------------------------------------------
// Convert real coordinates to integers.
//------------------------------------------------------------------------------
template<typename IntType, typename RealType>
void
QuantizedTessellation3d<IntType, RealType>::
quantize(const RealType* realcoords, IntType* intcoords) const {
  const RealType dx = length/(coordMax - coordMin);
  intcoords[0] = coordMin + IntType((realcoords[0] - xmin[0])/dx);
  intcoords[1] = coordMin + IntType((realcoords[1] - xmin[1])/dx);
  intcoords[2] = coordMin + IntType((realcoords[2] - xmin[2])/dx);
}

//------------------------------------------------------------------------------
// Convert int coordinates to reals.
//------------------------------------------------------------------------------
template<typename IntType, typename RealType>
void
QuantizedTessellation3d<IntType, RealType>::
dequantize(const IntType* intcoords, RealType* realcoords) const {
  const RealType dx = length/(coordMax - coordMin);
  realcoords[0] = xmin[0] + (intcoords[0] - coordMin)*dx;
  realcoords[1] = xmin[1] + (intcoords[1] - coordMin)*dx;
  realcoords[2] = xmin[2] + (intcoords[2] - coordMin)*dx;
}

//------------------------------------------------------------------------------
// Read out the current quantized tessellation state to a regular tessellation.
//------------------------------------------------------------------------------
template<typename IntType, typename RealType>
void
QuantizedTessellation3d<IntType, RealType>::
fillTessellation(Tessellation<3, RealType>& mesh) const {
  POLY_ASSERT2(false, "Implement me!");
}

//------------------------------------------------------------------------------
// Instantiate versions we know we need.
//------------------------------------------------------------------------------
template<typename IntType, typename RealType> 
IntType
QuantizedTessellation3d<IntType, RealType>::coordMin = std::numeric_limits<IntType>::min()/3;

template<typename IntType, typename RealType> 
IntType
QuantizedTessellation3d<IntType, RealType>::coordMax = std::numeric_limits<IntType>::max()/3;

template struct QuantizedTessellation3d<int, double>;
}
