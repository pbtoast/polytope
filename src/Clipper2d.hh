#ifndef POLYTOPE_CLIPPER2D_HH
#define POLYTOPE_CLIPPER2D_HH

#include <vector>

#include "Point.hh"

// We use the Boost.Geometry library to handle polygon intersections and such.
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/geometries/register/point.hpp>

namespace polytope {

template<typename RealType>
class Clipper2d {
public:

  // Typedefs for point, ring, and polygon type.
  // NOTE: Only Boost.Geometry exists for performing cell intersections.
  //       In the future, it may be necessary to template this class onto
  //       each of these types. After all, rings are just vectors of 
  //       points and polygons are just vectors of rings, specialized to
  //       outer points (facets) and inner points (holes). The ReducedPLC
  //       class -almost- handles everything we need.
  typedef Point2<RealType> PointType;
  typedef boost::geometry::model::ring   <PointType, false> RingType;
  typedef boost::geometry::model::polygon<PointType, false> PolygonType;

  // Constructor
  Clipper2d(PolygonType& boundary):
     mBoundary(boundary) {
  };

  // Destructor
  ~Clipper2d() {
  };

  // Boost.Geometry routine for performing the cell intersections. If more
  // than one ring is returned, we locate the one containing the initial 
  // generator and store any remainder pieces
  void clipCell(const PointType& cellGenerator,
                RingType& cellRing,
                std::vector<RingType>& orphans) {
    bool inside;
    unsigned j, k;
    std::vector<RingType> cellIntersections;
    boost::geometry::intersection(mBoundary, cellRing, cellIntersections);
    if (cellIntersections.size() == 0) {
       std::cerr << cellGenerator.x << " " << cellGenerator.y << std::endl 
                 << boost::geometry::dsv(cellRing) << std::endl
                 << boost::geometry::dsv(mBoundary) << std::endl;
    }
    POLY_ASSERT(cellIntersections.size() > 0);
    
    // Only one ring - this is the new celll
    if (cellIntersections.size() == 1) {
      cellRing = cellIntersections[0];
    } 
    // The ring containing the generator is the new cell
    else {
      k = cellIntersections.size();
      for (j = 0; j != cellIntersections.size(); ++j) {
        inside = boost::geometry::within(cellGenerator, cellIntersections[j]);
        if(inside)  k = j;
	else        orphans.push_back( cellIntersections[j] );
      }
      // If none of the cell intersections contain the point, check if
      // it's a vertex of the boundary
      if (k == cellIntersections.size()) {
        for (j = 0; j != cellIntersections.size(); ++j) {
          bool onBoundary = false;
          for (typename RingType::const_iterator itr = cellIntersections[j].begin();
               itr != cellIntersections[j].end(); ++itr) {
            onBoundary += (cellGenerator == *itr);
          }
          if(onBoundary) k = j;
        }
      }
      POLY_ASSERT(k < cellIntersections.size());
      cellRing = cellIntersections[k];
    }
  }
  
  // Store a reference to the boundary polygon
  PolygonType& mBoundary;

};

} //end namespace polytope

#endif
