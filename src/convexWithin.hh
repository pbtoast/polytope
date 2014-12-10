#ifndef POLYTOPE_CONVEXWITHIN_HH
#define POLYTOPE_CONVEXWITHIN_HH
//------------------------------------------------------------------------------
// convexWithin - determine if one convex polygon/polyhedron lies completely
// inside another one. We express each shape as a ReducedPLC for simplicity.
//------------------------------------------------------------------------------
#include <vector>
#include <limits>

#include "ReducedPLC.hh"
#include "polytope_geometric_utilities.hh"

namespace polytope {

//------------------------------------------------------------------------------
// Check if convex polygon a is within convex polygon b
//------------------------------------------------------------------------------
template<typename RealType>
bool
convexWithin(const ReducedPLC<2, RealType>& a, const ReducedPLC<2, RealType>& b) {
  const unsigned nvb = b.points.size() / 2;
  const unsigned nfb = b.facets.size();
  POLY_CONTRACT_VAR(nvb);

  bool inside = true;
  unsigned i, j, ifacet;

  ifacet = 0;
  while (inside and ifacet < nfb) {
    i = b.facets[ifacet][0];
    j = b.facets[ifacet][1];
    POLY_ASSERT(i < nvb);
    POLY_ASSERT(j < nvb);
    inside = (geometry::aboveBelow(b.points[2*i], b.points[2*i+1],
                                   b.points[2*j], b.points[2*j+1],
                                   a.points) == -1);
    ++ifacet;
  }
  return inside;
}

//------------------------------------------------------------------------------
// Check if convex polyhedron a is within convex polyhedron b
//------------------------------------------------------------------------------
template<typename RealType>
bool
convexWithin(const ReducedPLC<3, RealType>& a, const ReducedPLC<3, RealType>& b) {
  const unsigned nvb = b.points.size() / 2;
  const unsigned nfb = b.facets.size();
  POLY_CONTRACT_VAR(nvb);

  bool inside = true;
  unsigned i, j, k, ifacet;
  double nx, ny, nz;

  ifacet = 0;
  while (inside and ifacet < nfb) {
    i = b.facets[ifacet][0];
    j = b.facets[ifacet][1];
    k = b.facets[ifacet][2];
    POLY_ASSERT(i < nvb);
    POLY_ASSERT(j < nvb);
    POLY_ASSERT(k < nvb);
    geometry::computeNormal(b.points[3*i], b.points[3*i + 1], b.points[3*i + 2],
                            b.points[3*j], b.points[3*j + 1], b.points[3*j + 2],
                            b.points[3*k], b.points[3*k + 1], b.points[3*k + 2],
                            nx, ny, nz);
    inside = (geometry::aboveBelow(b.points[3*i], b.points[3*i+1], b.points[3*i+2],
                                   nx, ny, nz,
                                   a.points) == -1);
    ++ifacet;
  }
  return inside;
}

}

#endif
