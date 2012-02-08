//----------------------------------------------------------------------------//
// 2D implementation of the convex hull algorithm.
// Based on an example at http://www.algorithmist.com/index.php/Monotone_Chain_Convex_Hull.cpp
//----------------------------------------------------------------------------//
#ifndef __polytope_convexHull_2d__
#define __polytope_convexHull_2d__

#include <iostream>
#include <iterator>
#include <algorithm>

#include "polytope.hh"

namespace polytope {

//------------------------------------------------------------------------------
// A integer version of the simple 2D point.
//------------------------------------------------------------------------------
template<typename UintType>
struct Point2 {
  UintType x, y;
  Point2(): x(0), y(0) {}
  Point2(const UintType& xi, const UintType& yi): x(xi), y(yi) {}
  Point2& operator=(const Point2& rhs) { x = rhs.x; y = rhs.y; return *this; }
  bool operator==(const Point2& rhs) const { return (x == rhs.x and y == rhs.y); }
  bool operator<(const Point2& rhs) const {
    return (x < rhs.x                ? true :
            x == rhs.x and y < rhs.y ? true :
            false);
  }
  template<typename RealType>
  Point2(const RealType& xi, const RealType& yi, const RealType& dx): 
    x(static_cast<UintType>(xi/dx + 0.5)),
    y(static_cast<UintType>(yi/dx + 0.5)) {}
  template<typename RealType> RealType realx(const RealType& xmin, const RealType& dx) { return static_cast<RealType>(x*dx) + xmin; }
  template<typename RealType> RealType realy(const RealType& ymin, const RealType& dy) { return static_cast<RealType>(y*dy) + ymin; }
};

// It's nice being able to print these things.
template<typename UintType>
std::ostream&
operator<<(std::ostream& os, const Point2<UintType>& p) {
  os << "(" << p.x << " " << p.y << ")";
  return os;
}

namespace { // We hide internal functions in an anonymous namespace.

//------------------------------------------------------------------------------
// sign of the Z coordinate of cross product : (p2 - p1)x(p3 - p1).
//------------------------------------------------------------------------------
template<typename RealType>
int zcross_sign(const Point2<RealType>& p1, const Point2<RealType>& p2, const Point2<RealType>& p3) {
  double scale = 1.0/max(RealType(1), max(p1.x, max(p1.y, max(p2.x, max(p2.y, max(p3.x, p3.y))))));
  double ztest = (((p2.x - p1.x)*scale)*((p3.y - p1.y)*scale) -
                  ((p2.y - p1.y)*scale)*((p3.x - p1.x)*scale));
  return (ztest < 0.0 ? -1 :
          ztest > 0.0 ?  1 :
                         0);
  // return (p2.x - p1.x)*(p3.y - p1.y) - (p2.y - p1.y)*(p3.x - p1.x);
}

//------------------------------------------------------------------------------
// Comparator to compare std::pair's by their first element.
//------------------------------------------------------------------------------
template<typename T1, typename T2>
struct ComparePairByFirstElement {
  bool operator()(const std::pair<T1, T2>& lhs, const std::pair<T1, T2>& rhs) const {
    return lhs.first < rhs.first;
  }
};

} // end anonymous namespace

//------------------------------------------------------------------------------
// The method itself.
//------------------------------------------------------------------------------
template<typename RealType>
PLC<2, RealType>
convexHull_2d(const std::vector<RealType>& points,
              const RealType& xmin, const RealType& ymin,
              const RealType& dx) {
  typedef uint64_t CoordHash;
  typedef Point2<CoordHash> PointHash;

  ASSERT(points.size() % 2 == 0);
  const unsigned n = points.size() / 2;
  int i, j, k, t;
  
  // Hash the input points and sort them by x coordinate, remembering their original indices
  // in the input set.
  std::vector<std::pair<PointHash, unsigned> > sortedPoints;
  sortedPoints.reserve(n);
  for (i = 0; i != n; ++i) sortedPoints.push_back(std::make_pair(PointHash(CoordHash((points[2*i]     - xmin)/dx + 0.5),
                                                                           CoordHash((points[2*i + 1] - ymin)/dx + 0.5)),
                                                                 i));
  sort(sortedPoints.begin(), sortedPoints.end());

  // Prepare the result.
  std::vector<int> result(2*n);

  // Build the lower hull.
  for (i = 0, k = 0; i < n; i++) {
    while (k >= 2 and
           zcross_sign(sortedPoints[result[k - 2]].first, sortedPoints[result[k - 1]].first, sortedPoints[i].first) <= 0) k--;
    result[k++] = i;
  }

  // Build the upper hull.
  for (i = n - 2, t = k + 1; i >= 0; i--) {
    while (k >= t and
           zcross_sign(sortedPoints[result[k - 2]].first, sortedPoints[result[k - 1]].first, sortedPoints[i].first) <= 0) k--;
    result[k++] = i;
  }
  ASSERT(k >= 4);
  ASSERT(result.front() == result.back());

  // Translate our sorted information to a PLC based on the input point ordering and we're done.
  PLC<2, RealType> plc;
  for (i = 0; i != k - 1; ++i) {
    j = (i + 1) % k;
    plc.facets.push_back(std::vector<int>());
    plc.facets.back().push_back(sortedPoints[result[i]].second);
    plc.facets.back().push_back(sortedPoints[result[j]].second);
  }
  ASSERT(plc.facets.size() == k - 1);
  return plc;
}

}

#endif
