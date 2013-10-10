//----------------------------------------------------------------------------//
// 3D implementation of the convex hull algorithm.
// My own (possibly loose) interpretation of quick hull in 3D.
//----------------------------------------------------------------------------//
#ifndef __polytope_convexHull_3d__
#define __polytope_convexHull_3d__

#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <cmath>

#include "PLC.hh"
#include "polytope_internal.hh"
#include "polytope_geometric_utilities.hh"
#include "Point.hh"

namespace polytope {

namespace {

//------------------------------------------------------------------------------
// Distance from a point to a line.
// The line is specified by a point on the line (line0) and a unit direction
// vector (lineDirection).
//------------------------------------------------------------------------------
template<typename RealType>
RealType
pointLineDistance(const Point3<RealType>& p,
                  const Point3<RealType>& line0,
                  const Point3<RealType>& lineDirection) {
  const Point3<RealType> delta = p - line0;
  const RealType a = geometry::dot<3, RealType>(&delta.x, &lineDirection.x);
  const RealType c2 = geometry::dot<3, RealType>(&delta.x, &delta.x);
  return sqrt(std::abs(c2 - a*a));
}

//------------------------------------------------------------------------------
// Signed distance from a point to a plane
// The plane is specified by a point on the plane (line0) and a unit normal to
// the plane (planeNormal).
// The return value is signed based on if the point is above (> 0) or below (< 0)
// the plane based on the normal direction.
//------------------------------------------------------------------------------
template<typename RealType>
RealType
pointPlaneDistance(const Point3<RealType>& p,
                   const Point3<RealType>& plane0,
                   const Point3<RealType>& planeNormal) {
  const Point3<RealType> delta = p - plane0;
  return geometry::dot<3, RealType>(&delta.x, &planeNormal.x);
}

//------------------------------------------------------------------------------
// Comparison by (x,y, or z) coordinate.
//------------------------------------------------------------------------------
template<typename RealType>
bool compareMinX(const Point3<RealType>& a, const Point3<RealType>& b) {
  return a.x < b.x;
}

template<typename RealType>
bool compareMinY(const Point3<RealType>& a, const Point3<RealType>& b) {
  return a.y < b.y;
}

template<typename RealType>
bool compareMinZ(const Point3<RealType>& a, const Point3<RealType>& b) {
  return a.z < b.z;
}

//------------------------------------------------------------------------------
// Compare points by distance from a line.
//------------------------------------------------------------------------------
template<typename RealType>
struct CompareLineDistance {
  Point3<RealType> point, direction;
  CompareLineDistance(const Point3<RealType>& p,
                      const Point3<RealType>& dir): point(p), direction(dir) {}
  bool operator()(const Point3<RealType>& a, const Point3<RealType>& b) {
    const RealType adist = pointLineDistance(a, point, direction),
                   bdist = pointLineDistance(b, point, direction);
    return adist < bdist;
  }
};

//------------------------------------------------------------------------------
// Compare points by distance from a plane.
//------------------------------------------------------------------------------
template<typename RealType>
struct CompareAbsPlaneDistance {
  Point3<RealType> point, normal;
  CompareAbsPlaneDistance(const Point3<RealType>& p,
                          const Point3<RealType>& norm): point(p), normal(norm) {}
  bool operator()(const Point3<RealType>& a, const Point3<RealType>& b) {
    const RealType adist = std::abs(pointPlaneDistance(a, point, normal)),
                   bdist = std::abs(pointPlaneDistance(b, point, normal));
    return adist < bdist;
  }
};

//------------------------------------------------------------------------------
// Borrow the Point3 type as a tuple to create 3 node facets hashes.
//------------------------------------------------------------------------------
Point3<unsigned>
hashFacet(const unsigned i, const unsigned j, const unsigned k) {
  typedef Point3<unsigned> Tuple3;
  POLY_ASSERT(i != j and i != k and j != k);
  if (i < j and i < k) {
    if (j < k) {
      return Tuple3(i, j, k);
    } else {
      return Tuple3(i, k, j);
    }
  } else if (j < i and j < k) {
    if (i < k) {
      return Tuple3(j, i, k);
    } else {
      return Tuple3(j, k, i);
    }
  } else {
    if (i < j) {
      return Tuple3(k, i, j);
    } else {
      return Tuple3(k, j, i);
    }
  }
}

//------------------------------------------------------------------------------
// Hold one of our triangular facets.
//------------------------------------------------------------------------------
template<typename RealType>
struct TriangleFacet {
  typedef Point3<RealType> PointType;
  typedef Point3<unsigned> HashType;

  unsigned inode, jnode, knode;          // Indices of the vertices of the triangle.
  PointType normal;                      // Normal to the surface.
  std::vector<unsigned> remainingPoints; // Indices of points above this facet. 
  const std::vector<PointType>* points;  // The full set of point positions (i, j, k) refer to.

  TriangleFacet(): inode(0), jnode(0), knode(0), normal(), remainingPoints(), points(0) {}
  TriangleFacet(const unsigned i, 
                const unsigned j,
                const unsigned k,
                const std::vector<unsigned>& remainingPoints,
                const std::vector<PointType>& points): inode(i), 
                                                       jnode(j),
                                                       knode(k),
                                                       normal(),
                                                       remainingPoints(remainingPoints),
                                                       points(&points) {
    computeNormal();
  }

  void computeNormal() {
    const PointType ij = (*points)[jnode] - (*points)[inode];
    const PointType ik = (*points)[knode] - (*points)[inode];
    geometry::cross<3, RealType>(&ij.x, &ik.x, &normal.x);
    geometry::unitVector<3, RealType>(&normal.x);
  }

  int compare(const PointType& point) const { 
    const PointType ip = point - (*points)[inode];
    const RealType test = geometry::dot<3, RealType>(&ip.x, &normal.x);
    return (test < 0 ? -1 :
            test > 0 ?  1 :
            0);
  }

  void downselectRemainingPoints() {
    std::vector<unsigned> result;
    for (typename std::vector<unsigned>::const_iterator itr = remainingPoints.begin();
         itr != remainingPoints.end();
         ++itr) {
      if (this->compare((*points)[*itr]) == 1) result.push_back(*itr);
    }
    remainingPoints = result;
  }

  unsigned highestRemainingPoint() const {
    int result = -1;
    RealType maxAltitude = 0.0;
    for (typename std::vector<unsigned>::const_iterator itr = remainingPoints.begin();
         itr != remainingPoints.end();
         ++itr) {
      const RealType altitude = pointPlaneDistance((*points)[*itr], (*points)[inode], normal);
      if (altitude > maxAltitude) {
        maxAltitude = altitude;
        result = *itr;
      }
    }
    POLY_ASSERT(result >= 0);
    return unsigned(result);
  }

  void flip() { std::swap(inode, jnode); normal *= -1; }
  bool operator==(const TriangleFacet& rhs) const { return hashFacet(inode, jnode, knode) == hashFacet(rhs.inode, rhs.jnode, rhs.knode); }
  bool operator<(const TriangleFacet& rhs) const { return hashFacet(inode, jnode, knode) < hashFacet(rhs.inode, rhs.jnode, rhs.knode); }

  // output operator
  friend std::ostream& operator<<(std::ostream& os, const TriangleFacet& facet) {
    os << "TriangleFacet(" << facet.inode << " " << facet.jnode << " " << facet.knode 
       << ") : remaining points = (";
    for (std::vector<unsigned>::const_iterator itr = facet.remainingPoints.begin();
         itr != facet.remainingPoints.end();
         ++itr) os << " " << *itr;
    os << ")";
    return os;
  }
};

}

//------------------------------------------------------------------------------
// The 3D convex hull itself.  This is the one users should call -- it forwards
// all work to the worker ConvexHull3d class.
//------------------------------------------------------------------------------
template<typename RealType>
PLC<3, RealType>
convexHull_3d(const std::vector<RealType>& points,
              const RealType* low,
              const RealType& dx) {
  typedef Point3<RealType> PointType;

  // Pre-conditions.
  POLY_ASSERT(points.size() % 3 == 0);

  // Reduce to the unique set of input points.
  std::vector<RealType> upoints;
  std::vector<unsigned> pointMap;
  geometry::uniquePoints<3, RealType>(points, upoints, pointMap);
  const unsigned nunique = upoints.size()/3;
  POLY_ASSERT(nunique >= 4);
  std::vector<PointType> ps(nunique);
  for (unsigned i = 0; i != nunique; ++i) {
    ps[i].x = upoints[3*i];
    ps[i].y = upoints[3*i+1];
    ps[i].z = upoints[3*i+2];
    ps[i].index = i;
  }

  // Find the points on the min/max x, y, & z.  These must be in the hull.
  const unsigned ixmin = std::min_element(ps.begin(), ps.end(), compareMinX<RealType>)->index,
                 ixmax = std::max_element(ps.begin(), ps.end(), compareMinX<RealType>)->index,
                 iymin = std::min_element(ps.begin(), ps.end(), compareMinY<RealType>)->index,
                 iymax = std::max_element(ps.begin(), ps.end(), compareMinY<RealType>)->index,
                 izmin = std::min_element(ps.begin(), ps.end(), compareMinZ<RealType>)->index,
                 izmax = std::max_element(ps.begin(), ps.end(), compareMinZ<RealType>)->index;
  POLY_ASSERT(ixmin != ixmax and
              iymin != iymax and
              izmin != izmax);

  // Choose the two points furthest apart.
  const PointType boxx = ps[ixmax] - ps[ixmin], boxy = ps[iymax] - ps[iymin], boxz = ps[izmax] - ps[izmin];
  const RealType boxx2 = geometry::dot<3, RealType>(&boxx.x, &boxx.x), 
                 boxy2 = geometry::dot<3, RealType>(&boxy.x, &boxy.x),
                 boxz2 = geometry::dot<3, RealType>(&boxz.x, &boxz.x);
  POLY_ASSERT(boxx2 > 0.0 and boxy2 > 0.0 and boxz2 > 0.0);
  TriangleFacet<RealType> startingFacet;
  startingFacet.points = &ps;
  if (boxx2 >= std::max(boxy2, boxz2)) {
    startingFacet.inode = ixmin;
    startingFacet.jnode = ixmax;
  } else if (boxy2 >= boxz2) {
    startingFacet.inode = iymin;
    startingFacet.jnode = iymax;
  } else {
    POLY_ASSERT(boxz2 >= std::max(boxx2, boxy2));
    startingFacet.inode = izmin;
    startingFacet.jnode = izmax;
  }

  // Select a third point as the most distant from the line defined by the two we just picked.
  // This one should also be on the hull.
  PointType lineDirection = ps[startingFacet.jnode] - ps[startingFacet.inode];
  geometry::unitVector<3, RealType>(&lineDirection.x);
  startingFacet.knode = std::max_element(ps.begin(), ps.end(), CompareLineDistance<RealType>(ps[startingFacet.inode], lineDirection))->index;
  startingFacet.computeNormal();
  startingFacet.remainingPoints.resize(ps.size());
  for (unsigned i = 0; i != ps.size(); ++i) startingFacet.remainingPoints[i] = i;

  // We have the triangular base, so pick the furthest point from this base as the apex of our
  // starting tetrahedron.  This point should also be on the hull.
  const unsigned apex = std::max_element(ps.begin(), ps.end(), CompareAbsPlaneDistance<RealType>(ps[startingFacet.inode], startingFacet.normal))->index;

  // If the apex point is below our starting facet, flip the facet.
  if (startingFacet.compare(ps[apex]) == 1) startingFacet.flip();

  // Create our working facets from the starting tetrahedron.
  typedef typename std::vector<TriangleFacet<RealType> > FacetSet;
  FacetSet workingFacets;
  workingFacets.push_back(startingFacet);
  workingFacets.push_back(TriangleFacet<RealType>(startingFacet.jnode, startingFacet.inode, apex, startingFacet.remainingPoints, ps));
  workingFacets.push_back(TriangleFacet<RealType>(startingFacet.inode, startingFacet.knode, apex, startingFacet.remainingPoints, ps));
  workingFacets.push_back(TriangleFacet<RealType>(startingFacet.knode, startingFacet.jnode, apex, startingFacet.remainingPoints, ps));
  for (typename FacetSet::iterator facetItr = workingFacets.begin();
       facetItr != workingFacets.end();
       ++facetItr) {
    facetItr->downselectRemainingPoints();
  }

  // Iterate until all points are interior to our facet set.
  bool done = false;
  while (not done) {
    done = true;
    FacetSet newFacets;

    std::cerr << "NEW PASS!" << std::endl;
    for (typename FacetSet::iterator facetItr = workingFacets.begin();
         facetItr != workingFacets.end();
         ++facetItr) {

      // Is this facet done yet?
      if (facetItr->remainingPoints.empty()) {

        // Yep, copy it to the new set.
        std::cerr << "      Final facet : " << (*facetItr) << std::endl;
        newFacets.push_back(*facetItr);

      } else {

        // Nope, we need to refine new facets from this one.
        done = false;

        // Find the point furthest above this facet.
        const unsigned apex = facetItr->highestRemainingPoint();
        std::cerr << "      Refining facet : " << (*facetItr) << " with new apex " << apex << std::endl;

        // Add our new facets.
        TriangleFacet<RealType> afacet(startingFacet.jnode, startingFacet.inode, apex, facetItr->remainingPoints, ps),
                                bfacet(startingFacet.inode, startingFacet.knode, apex, facetItr->remainingPoints, ps),
                                cfacet(startingFacet.knode, startingFacet.jnode, apex, facetItr->remainingPoints, ps);
        afacet.downselectRemainingPoints();
        bfacet.downselectRemainingPoints();
        cfacet.downselectRemainingPoints();
        std::cerr << "   new facet : " << afacet << std::endl
                  << "   new facet : " << bfacet << std::endl
                  << "   new facet : " << cfacet << std::endl;
        newFacets.push_back(afacet);
        newFacets.push_back(bfacet);
        newFacets.push_back(cfacet);
      }
    }

    // Assign the new set of Facets.
    workingFacets = newFacets;
  }

  // Read out the data to the PLC and we're done.
  PLC<3, RealType> plc;
  for (typename FacetSet::iterator facetItr = workingFacets.begin();
       facetItr != workingFacets.end();
       ++facetItr) {
    plc.facets.push_back(std::vector<int>());
    plc.facets.back().push_back(pointMap[facetItr->inode]);
    plc.facets.back().push_back(pointMap[facetItr->jnode]);
    plc.facets.back().push_back(pointMap[facetItr->knode]);
  }
  return plc;
}

}

#endif
