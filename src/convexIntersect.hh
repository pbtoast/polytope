#ifndef __polytope_convexIntersection__
#define __polytope_convexIntersection__

#include <vector>

#include "PLC.hh"

namespace { // anonymous

//------------------------------------------------------------------------------
// Compare a point to a line segment, and determine if the point is above,
// below, or colinear with the line (assuming the line points (l1, l2) are
// specified in counter-clockwise manner to indicate the interior direction.
//------------------------------------------------------------------------------
template<typename RealType>
int compare(const RealType& l1x, const RealType& l1y,
            const RealType& l2x, const RealType& l2y,
            const RealType& px, const RealType& py) {
  // We scale here to avoid potential overflow issues with integer types.
  double scale = 1.0/max(RealType(1), max(l1x, max(l1y, max(l2x, max(l2y, max(px, py))))));
  double ztest = (((l2x - l1x)*scale)*((py - l1y)*scale) -
                  ((l2y - l1y)*scale)*((px - l1x)*scale));
  return (ztest < 0.0 ? -1 :
          ztest > 0.0 ?  1 :
                         0);
}

}           // anonymous

namespace polytope {

//------------------------------------------------------------------------------
// Convex polygon intersection.
//------------------------------------------------------------------------------
template<typename RealType>
bool
convexIntersect(const PLC<2, RealType>& a, const PLC<2, RealType>& b,
                const std::vector<RealType>& apoints, const std::vector<RealType>& bpoints) {
  const unsigned nva = apoints.size() / 2;
  const unsigned nvb = bpoints.size() / 2;

  unsigned i, k, ia, ja, ib, jb;

  // Check if we can exclude b from a.
  bool outside = true;
  {
    i = 0;
    while (outside and i < nva) {
      ia = a.facets[i][0];
      ja = a.facets[i][1];
      k = 0;
      while (outside and k != nvb) {
        ib = b.facets[k][0];
        outside = (compare(apoints[2*ia], apoints[2*ia + 1],
                           apoints[2*ja], apoints[2*ja + 1],
                           bpoints[2*ib], bpoints[2*ib + 1]) == 1);
        ++k;
      }
      ++i;
    }
    if (outside) return false;
  }

  // Check if we can exclude a from b.
  outside = true;
  {
    i = 0;
    while (outside and i < nvb) {
      ib = b.facets[i][0];
      jb = b.facets[i][1];
      k = 0;
      while (outside and k != nva) {
        ia = a.facets[k][0];
        outside = (compare(bpoints[2*ib], bpoints[2*ib + 1],
                           bpoints[2*jb], bpoints[2*jb + 1],
                           apoints[2*ia], apoints[2*ia + 1]) == 1);
        ++k;
      }
      ++i;
    }
    if (outside) return false;
  }

  // We can't exclude anybody, so must intersect!
  return true;
}

}

#endif
