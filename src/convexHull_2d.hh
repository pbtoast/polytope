//----------------------------------------------------------------------------//
// 2D implementation of the convex hull algorithm.
// Based on an example at http://www.algorithmist.com/index.php/Monotone_Chain_Convex_Hull.cpp
//----------------------------------------------------------------------------//
#ifndef __polytope_convexHull_2d__
#define __polytope_convexHull_2d__

#include <iostream>
#include <iterator>
#include <algorithm>
#include <set>
#include <stdint.h>

#include "polytope.hh"
#include "Point.hh"

namespace polytope {

namespace { // We hide internal functions in an anonymous namespace.

//------------------------------------------------------------------------------
// A fuzzy comparison operator for our quantized Point2 type.
//------------------------------------------------------------------------------
template<typename UintType>
struct FuzzyPoint2LessThan {
  UintType fuzz;
  FuzzyPoint2LessThan(const UintType ifuzz = 1): fuzz(ifuzz) {}
  bool operator()(const Point2<UintType>& p1, const Point2<UintType>& p2) {
    return (int(p2.x) - int(p1.x) > fuzz ? true :
            int(p2.y) - int(p1.y) > fuzz ? true : false);
  }
  bool operator()(const std::pair<Point2<UintType>, unsigned>& p1,
                  const std::pair<Point2<UintType>, unsigned>& p2) {
    return operator()(p1.first, p2.first);
  }
};

//------------------------------------------------------------------------------
// sign of the Z coordinate of cross product : (p2 - p1)x(p3 - p1).
//------------------------------------------------------------------------------
template<typename RealType>
int zcross_sign(const Point2<RealType>& p1, const Point2<RealType>& p2, const Point2<RealType>& p3) {
//   double scale = 1.0/max(RealType(1), max(p1.x, max(p1.y, max(p2.x, max(p2.y, max(p3.x, p3.y))))));
  const double ztest = 
    (double(p2.x) - double(p1.x))*(double(p3.y) - double(p1.y)) -
    (double(p2.y) - double(p1.y))*(double(p3.x) - double(p1.x));
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
              const RealType* low,
              const RealType& dx) {
  typedef uint64_t CoordHash;
  typedef Point2<CoordHash> PointHash;

  ASSERT(points.size() % 2 == 0);
  const unsigned n = points.size() / 2;
  int i, j, k, t;

  // Hash the input points and sort them by x coordinate, remembering their original indices
  // in the input set.  We also ensure that only unique (using a fuzzy comparison) points
  // are inserted here, since duplicates mess up the hull calculation.
  const RealType& xmin = low[0];
  const RealType& ymin = low[1];
  std::set<std::pair<PointHash, unsigned>, FuzzyPoint2LessThan<CoordHash> > uniquePoints;
  for (i = 0; i != n; ++i) {
    uniquePoints.insert(std::make_pair(PointHash(CoordHash((points[2*i]     - xmin)/dx + 0.5),
                                                 CoordHash((points[2*i + 1] - ymin)/dx + 0.5)),
                                       i));
  }
  std::vector<std::pair<PointHash, unsigned> > sortedPoints(uniquePoints.begin(), uniquePoints.end());
  std::sort(sortedPoints.begin(), sortedPoints.end());

//   std::vector<std::pair<PointHash, unsigned> > sortedPoints;
//   sortedPoints.reserve(n);
//   for (i = 0; i != n; ++i) {
//     PointHash p(CoordHash((points[2*i]     - xmin)/dx + 0.5),
//                 CoordHash((points[2*i + 1] - ymin)/dx + 0.5));
//     bool uniquePoint = true;
//     j = 0;
//     while (uniquePoint and j != sortedPoints.size()) {
//       uniquePoint = (abs(p.x - sortedPoints[j].first.x) > 1 or
//                      abs(p.y - sortedPoints[j].first.y) > 1);
//       ++j;
//     }
//     if (uniquePoint) sortedPoints.push_back(std::make_pair(p, i));
//   }

  // Prepare the result.
  const unsigned nunique = sortedPoints.size();
  std::vector<int> result(2*nunique);

  // Build the lower hull.
  for (i = 0, k = 0; i < nunique; i++) {
    while (k >= 2 and
           zcross_sign(sortedPoints[result[k - 2]].first, sortedPoints[result[k - 1]].first, sortedPoints[i].first) <= 0) k--;
    result[k++] = i;
  }

  // Build the upper hull.
  for (i = nunique - 2, t = k + 1; i >= 0; i--) {
    while (k >= t and
           zcross_sign(sortedPoints[result[k - 2]].first, sortedPoints[result[k - 1]].first, sortedPoints[i].first) <= 0) k--;
    result[k++] = i;
  }
  if (!(k >= 4)) {
    std::cerr << "Blago!  " << n << " " << nunique << " " << k << std::endl;
    for (unsigned i = 0; i != nunique; ++i) std::cerr << "  --> " << sortedPoints[i].first << std::endl;
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
