//------------------------------------------------------------------------
// QuantizedTessellation2d
//
// An intermediate representation for 2D tessellations in integer
// coordinates.
//------------------------------------------------------------------------
#ifndef __Polytope_QuantizedTessellation2d__
#define __Polytope_QuantizedTessellation2d__

#include <vector>
#include <utility>   // For std::pair

#include "Point.hh"
#include "polytope_tessellator_utilities.hh"

namespace polytope {

template<typename IntType, typename RealType>
struct QuantizedTessellation2d {
  typedef Point2<IntType> IntPoint;
  typedef Point2<RealType> RealPoint;
  
  RealType xmin[2], xmax[2], length, infRadius;
  std::vector<IntPoint> generators;

  std::vector<IntPoint> nodes;
  std::vector<std::pair<int, int> > edges;
  std::vector<std::vector<int> > cellEdges;

  // Construct with the given generators.  Finds the bounding limits, sets the infRadius,
  // and sets the quantized generators.
  QuantizedTessellation2d(const std::vector<RealType>& points);

  // Convert real coordinates to integers.
  void quantize(const RealType* realcoords, IntType* intcoords) const;

  // Convert int coordinates to reals.
  void dequantize(const IntType* intcoords, RealType* realcoords) const;

  // Read out the current QuantizedTessellation to regular Tessellation.
  void fillTessellation(Tessellation<2, RealType>& mesh) const;

  // Static properties about our coordinates.
  static IntType coordMin, coordMax;
};

}

#endif
