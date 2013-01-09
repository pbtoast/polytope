#ifndef POLYTOPE_NEAREST_POINT_HH
#define POLYTOPE_NEAREST_POINT_HH
//------------------------------------------------------------------------------
// nearestPoint - compute the closest point on a PLC to a given point.
//
// The point can be on the interior or exterior of the PLC.  Checking of holes
// is also included.
// The point is returned as the last argument, and the distance as the result.
//------------------------------------------------------------------------------
#include <vector>
#include <limits>

#include "PLC.hh"
#include "polytope_internal.hh"

namespace polytope {

//------------------------------------------------------------------------------
// Local utility methods.
//------------------------------------------------------------------------------
namespace {

// Distance between points.
template<int Dimension, typename RealType> RealType distance(const RealType* a, const RealType* b);
template<typename RealType> RealType distance<2, RealType>(const RealType* a, const RealType* b) {
  const RealType dx = a[0] - b[0];
  const RealType dy = a[1] - b[1];
  return std::sqrt(dx*dx + dy*dy);
}
template<typename RealType> RealType distance<3, RealType>(const RealType* a, const RealType* b) {
  const RealType dx = a[0] - b[0];
  const RealType dy = a[1] - b[1];
  const RealType dz = a[2] - b[2];
  return std::sqrt(dx*dx + dy*dy + dz*dz);
}

// Dot product.
template<int Dimension, typename RealType> RealType dot(const RealType* a, const RealType* b);
template<typename RealType> RealType dot<2, RealType>(const RealType* a, const RealType* b) {
  return a[0]*b[0] + a[1]+b[1];
}
template<typename RealType> RealType dot<3, RealType>(const RealType* a, const RealType* b) {
  return a[0]*b[0] + a[1]+b[1] + a[2]*b[2];
}

// Find the closest point on a line segment (2D).
template<typename RealType> 
void
closestPointOnSegment2D(const RealType* point, 
                        const RealType* s1,
                        const RealType* s2,
                        RealType* result) {
  const RealType shat[2] = {s2[0] - s1[0], s2[1] - s1[1]};
  const RealType seglength = std::sqrt(shat[0]*shat[0] + shat[1]*shat[1]);
  if (seglength < 1.0e-10) {
    result[0] = 0.5*(s1[0] + s2[0]);
    result[1] = 0.5*(s1[1] + s2[1]);
  } else {
    shat[0] /= seglength;
    shat[1] /= seglength;
    const RealType s1p[2] = {point[0] - s1[0], point[1] - s1[1]};
    const RealType ptest = std::max(0.0, std::min(seglength, dot<2, RealType>(s1p, shat)));
    result[0] = s1[0] + ptest*shat[0];
    result[1] = s1[1] + ptest*shat[1];
  }
}

// Find the closest point on a given set of facets.
template<int Dimension, typename RealType> 
RealType
closestPointOnFacets(const RealType* point,
                     const unsigned numVertices,
                     const RealType* vertices,
                     const std::vector<std::vector<int> >& facets,
                     RealType* result);

// 2-D specialization.
template<typename RealType>
RealType
closestPointOnFacets<2, RealType>(const RealType* point,
                                  const unsigned numVertices,
                                  const RealType* vertices,
                                  const std::vector<std::vector<int> >& facets,
                                  RealType* result) {
  unsigned i, j;
  RealType dist, minDist = std::numeric_limits<RealType>::max();
  RealType[2] candidate;
  const unsigned numFacets = facets.size();
  for (unsigned ifacet = 0; ifacet != numFacets; ++ifacet) {
    POLY_ASSERT(facets[ifacet].size() == 2);
    i = facets[ifacet][0];
    j = facets[ifacet][1];
    POLY_ASSERT(i >= 0 and i < numVertices);
    POLY_ASSERT(j >= 0 and j < numVertices);
    closestPointOnSegment2D(point, &vertices[2*i], &vertices[2*j], candidate);
    dist = distance<2, RealType>(point, candidate);
    if (dist < minDist) {
      minDist = dist;
      result[0] = candidate[0];
      result[1] = candidate[1];
    }
  }
  return minDist;
}

}

//------------------------------------------------------------------------------
// nearestPoint
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
RealType
nearestPoint(const RealType* point,
             const unsigned numVertices,
             const RealType* vertices,
             const PLC<Dimension, RealType>& plc,
             RealType* result) {

  // Check the outer boundary of the PLC.
  RealType minDist = closestPointOnFacets<Dimension, RealType>(point, numVertices, vertices, plc.facets, result);

  // Check each of the holes.
  RealType dist, candidate[2];
  for (unsigned ihole = 0; ihole != plc.holes.size(); ++ihole) {
    dist = closestPointOnFacets<Dimension, RealType>(point, numVertices, vertices, plc.holes[ihole], candidate);
    if (dist < minDist) {
      minDist = dist;
      for (unsigned i = 0; i != Dimension; ++i) result[i] = candidate[i];
    }
  }
  return minDist;
}

}

#endif
