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

namespace polytope {

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

  typedef int64_t CoordHash;
  typedef convexHull_helpers::ConvexHull3d<CoordHash>::Point Point;

  // Pre-conditions.
  POLY_ASSERT(points.size() % 3 == 0);
  const unsigned n = points.size() / 3;

  // Reduce to the unique set of input points.
  std::vector<RealType> upoints;
  std::vector<unsigned> pointMap;
  geometry::uniquePoints<3, RealType>(points, uniquePoints, pointMap);
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
  const unsigned ixmin = std::min_element(ps.begin(), ps.end(), CompareMinX<RealType>())->index,
                 ixmax = std::max_element(ps.begin(), ps.end(), CompareMaxX<RealType>())->index,
                 iymin = std::min_element(ps.begin(), ps.end(), CompareMinY<RealType>())->index,
                 iymax = std::max_element(ps.begin(), ps.end(), CompareMaxY<RealType>())->index,
                 izmin = std::min_element(ps.begin(), ps.end(), CompareMinZ<RealType>())->index,
                 izmax = std::max_element(ps.begin(), ps.end(), CompareMaxZ<RealType>())->index;
  POLY_ASSERT(ixmin != ixmax and
              iymin != iymax and
              izmin != izmax);

  

  const RealType& xmin = low[0];
  const RealType& ymin = low[1];
  const RealType& zmin = low[2];

  // Convert the input coordinates to unique integer point types.  Simultaneously we 
  // reduce to the unique set of points.
  typedef std::set<Point, convexHull_helpers::ConvexHull3d<CoordHash>::FuzzyPointLessThan> Set;
  Set pointSet;
  for (i = 0; i != n; ++i) {
    pointSet.insert(Point(CoordHash((points[3*i]     - xmin)/dx + 0.5),
                          CoordHash((points[3*i + 1] - ymin)/dx + 0.5),
                          CoordHash((points[3*i + 2] - zmin)/dx + 0.5), 
                          i));
  }
  POLY_ASSERT(pointSet.size() <= n);

  // Extract the unique set of points to a vector.
  std::vector<Point> uniquePoints(pointSet.begin(), pointSet.end());
  POLY_ASSERT(uniquePoints.size() == pointSet.size());

  // // Blago!
  // std::cout << "Unique points: " << std::endl;
  // for (unsigned k = 0; k != uniquePoints.size(); ++k) {
  //   std::cout << "   ---> " << uniquePoints[k].x << " "<< uniquePoints[k].y << " " << uniquePoints[k].z << " " << uniquePoints[k].index << std::endl;
  // }
  // // Blago!

  // Dispatch the work 
  const std::vector<std::vector<unsigned> > faces = convexHull_helpers::ConvexHull3d<CoordHash>::process(uniquePoints);

  // Read out the data to the PLC and we're done.
  PLC<3, RealType> plc;
  // RealType xc, yc, zc;
  for (i = 0; i != faces.size(); ++i) {
    POLY_ASSERT(faces[i].size() == 3);
    plc.facets.push_back(std::vector<int>());
    plc.facets.back().push_back(faces[i][0]);
    plc.facets.back().push_back(faces[i][1]);
    plc.facets.back().push_back(faces[i][2]);
    // convexHull_helpers::faceCentroid(points, faces[i], xc, yc, zc);
    // std::cerr << "  -----> " << faces[i][0] << " " << faces[i][1] << " " << faces[i][2] << " : (" << xc << " " << yc << " " << zc << ")" << std::endl;
  }
  return plc;
}

}

#endif
