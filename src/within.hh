#ifndef POLYTOPE_WITHIN_HH
#define POLYTOPE_WITHIN_HH
//------------------------------------------------------------------------------
// within - determine whether a point lies inside a complex boundary
//
// Checks if a point is inside a PLC boundary AND outside holes, if present.
//------------------------------------------------------------------------------
#include <vector>
#include <limits>

#include "PLC.hh"
#include "polytope_internal.hh"
#include "polytope_geometric_utilities.hh"

namespace polytope {

//------------------------------------------------------------------------------
// Local utility methods.
//------------------------------------------------------------------------------
namespace {

// Determine if point lies inside a set of facets
// Functor definition first.
template<int Dimension, typename RealType> struct WithinBoundaryFunctor;

// 2-D specialization.
template<typename RealType>
struct WithinBoundaryFunctor<2, RealType> {
  static bool impl(const RealType* point,
                   const unsigned numVertices,
                   const RealType* vertices,
                   const std::vector<std::vector<int> >& facets) {
    const unsigned numFacets = facets.size();
    unsigned i = facets[0][0];
    bool isInside = geometry::withinPolygon2D( point, numFacets, &vertices[2*i] );
    return isInside;
  }
};

// Functional interface.
template<int Dimension, typename RealType> 
bool withinBoundary(const RealType* point,
                    const unsigned numVertices,
                    const RealType* vertices,
                    const std::vector<std::vector<int> >& facets) {
  return WithinBoundaryFunctor<2, RealType>::impl(point, numVertices, vertices, facets);
}

}

//------------------------------------------------------------------------------
// within
//------------------------------------------------------------------------------
template<int Dimension, typename RealType>
bool
within(const RealType* point,
       const unsigned numVertices,
       const RealType* vertices,
       const PLC<Dimension, RealType>& plc) {

  // Check the outer boundary of the PLC.
  bool isInside = withinBoundary<Dimension, RealType>(point, numVertices, vertices, plc.facets);

  // Check each of the holes.
  for (unsigned ihole = 0; ihole != plc.holes.size(); ++ihole) {
    isInside ^= withinBoundary<Dimension, RealType>(point, numVertices, vertices, plc.holes[ihole]);
  }
  return isInside;
}

}

#endif
