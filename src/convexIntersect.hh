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
  using std::max;
  using std::abs;
  double scale = 1.0/max(RealType(1), 
                         max(abs(l1x), max(abs(l1y),
                                           max(abs(l2x), max(abs(l2y),
                                                             max(abs(px), abs(py)))))));
  double ztest = (((l2x - l1x)*scale)*((py - l1y)*scale) -
                  ((l2y - l1y)*scale)*((px - l1x)*scale));
  return (ztest < 0.0 ? -1 :
          ztest > 0.0 ?  1 :
                         0);
}

//------------------------------------------------------------------------------
// Compare a point to a plane, and determine if the point is above,
// below, or coplanar with the plane.  The plane is specified by the 
// the point normal [(ox, oy, oz), (nx, ny, nz)], and the point by (px, py, pz).
//------------------------------------------------------------------------------
template<typename RealType>
int compare(const RealType& ox, const RealType& oy, const RealType& oz,
            const double& nx, const double& ny, const double& nz,
            const RealType& px, const RealType& py, const RealType& pz) {
  const double ztest = nx*double(px - ox) + ny*double(py - oy) + nz*double(pz - oz);
  return (ztest < 0.0 ? -1 :
          ztest > 0.0 ?  1 :
                         0);
}

//------------------------------------------------------------------------------
// Compute the 3D normal given three points (a, b, c).
//------------------------------------------------------------------------------
template<typename RealType>
void computeNormal(const RealType& ax, const RealType& ay, const RealType& az,
                   const RealType& bx, const RealType& by, const RealType& bz,
                   const RealType& cx, const RealType& cy, const RealType& cz,
                   double& nx, double& ny, double& nz) {
  const double dx_ab = bx - ax, dy_ab = by - ay, dz_ab = bz - az;
  const double dx_ac = cx - ax, dy_ac = cy - ay, dz_ac = cz - az;
  nx = dy_ab*dz_ac - dz_ab*dz_ac;
  ny = dz_ab*dx_ac - dx_ab*dz_ac;
  nz = dx_ab*dy_ac - dy_ab*dx_ac;
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

//------------------------------------------------------------------------------
// Convex polyhedron intersection.
//------------------------------------------------------------------------------
template<typename RealType>
bool
convexIntersect(const PLC<3, RealType>& a, const PLC<3, RealType>& b,
                const std::vector<RealType>& apoints, const std::vector<RealType>& bpoints) {
  const unsigned nva = apoints.size() / 3;
  const unsigned nvb = bpoints.size() / 3;
  const unsigned naf = a.facets.size();
  const unsigned nbf = b.facets.size();

  unsigned ifaceta, ifacetb, i, j, k;
  double nx, ny, nz;

  // Check if we can exclude b from a.
  bool outside = true;
  {
    ifaceta = 0;
    while (outside and ifaceta != naf) {
      ASSERT(a.facets[ifaceta].size() >= 3);
      i = a.facets[ifaceta][0];
      j = a.facets[ifaceta][1];
      k = a.facets[ifaceta][2];
      ASSERT(i < nva and j < nva and k < nva);
      computeNormal(apoints[3*i], apoints[3*i + 1], apoints[3*i + 2],
                    apoints[3*j], apoints[3*j + 1], apoints[3*j + 2],
                    apoints[3*k], apoints[3*k + 1], apoints[3*k + 2],
                    nx, ny, nz);
      ifacetb = 0;
      while (outside and ifacetb != nbf) {
        j = 0;
        while (outside and j != b.facets[ifacetb].size()) {
          k = b.facets[ifacetb][j];
          ASSERT(k < nvb);
          outside = (compare(apoints[3*i], apoints[3*i + 1], apoints[3*i + 2],
                             nx, ny, nz,
                             bpoints[3*k], bpoints[3*k + 1], bpoints[3*k + 2]) == 1);
          ++j;
        }
        ++ifacetb;
      }
      ++ifaceta;
    }
    if (outside) return false;
  }

  // Check if we can exclude a from b.
  {
    ifacetb = 0;
    while (outside and ifacetb != nbf) {
      ASSERT(b.facets[ifacetb].size() >= 3);
      i = b.facets[ifacetb][0];
      j = b.facets[ifacetb][1];
      k = b.facets[ifacetb][2];
      ASSERT(i < nvb and j < nvb and k < nvb);
      computeNormal(bpoints[3*i], bpoints[3*i + 1], bpoints[3*i + 2],
                    bpoints[3*j], bpoints[3*j + 1], bpoints[3*j + 2],
                    bpoints[3*k], bpoints[3*k + 1], bpoints[3*k + 2],
                    nx, ny, nz);
      ifaceta = 0;
      while (outside and ifaceta != naf) {
        j = 0;
        while (outside and j != a.facets[ifaceta].size()) {
          k = a.facets[ifaceta][j];
          ASSERT(k < nva);
          outside = (compare(bpoints[3*i], bpoints[3*i + 1], bpoints[3*i + 2],
                             nx, ny, nz,
                             apoints[3*k], apoints[3*k + 1], apoints[3*k + 2]) == 1);
          ++j;
        }
        ++ifaceta;
      }
      ++ifacetb;
    }
    if (outside) return false;
  }

  // We can't exclude anybody, so must intersect!
  return true;
}

}

#endif
