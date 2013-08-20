#ifndef POLYTOPE_CLIPPER2D_HH
#define POLYTOPE_CLIPPER2D_HH

#if HAVE_BOOST

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
     // Pre-condition
     POLY_ASSERT(!boost::geometry::intersects(boundary));
  }

  // Destructor
  ~Clipper2d() {}

  // Boost.Geometry routine for performing the cell intersections. If more
  // than one ring is returned, we locate the one containing the initial 
  // generator and store any remainder pieces
  void clipCell(const PointType& cellGenerator,
                RingType& cellRing,
                std::vector<RingType>& orphans) {

    bool inside, onBoundary;
    unsigned j, k;
    if (boost::geometry::intersects(cellRing, mBoundary)) {
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
	  // inside = boost::geometry::within(cellGenerator, cellIntersections[j]);
	  inside = boost::geometry::covered_by(cellGenerator, cellIntersections[j]);
	  if (inside) k = j;
	  else {
	    onBoundary = false;
	    for (typename RingType::const_iterator itr = cellIntersections[j].begin();
		 itr != cellIntersections[j].end(); 
		 ++itr)  onBoundary += (cellGenerator == *itr);
	    if (onBoundary) k = j;
	    else            orphans.push_back(cellIntersections[j]);
	  }
	}
	POLY_ASSERT(k < cellIntersections.size());
	cellRing = cellIntersections[k];
      }
     

//       // Flag any post-clipped points that coincide with the boundary
//       for (typename RingType::iterator itr = cellRing.begin();
// 	   itr != cellRing.end()-1; ++itr) {
// 	for (typename RingType::const_iterator oItr = mBoundary.outer().begin();
// 	     oItr != mBoundary.outer().end();
// 	     ++oItr) {
// 	  if (*itr == *oItr)  itr->index = 1;
// 	}
// 	typename std::vector<RingType>& holes = mBoundary.inners();
// 	for (unsigned ihole = 0; ihole != holes.size(); ++ihole) {
// 	  for (typename RingType::const_iterator iItr = holes[ihole].begin();
// 	       iItr != holes[ihole].end()-1;
// 	       ++iItr) {
// 	    if (*itr == *iItr)  itr->index = 1;
// 	  }
// 	}
	
// // 	//Blago!
// // 	if (itr->index == 1) std::cerr << (*itr) << std::endl;
// // 	//Blago!
//       }

 
    }

    // Post-conditions
    POLY_ASSERT(cellRing.size() > 2);
    POLY_ASSERT(cellRing.front() == cellRing.back());
    //POLY_ASSERT(!boost::geometry::intersects(cellRing));
  }
  
  // Store a reference to the boundary polygon
  PolygonType& mBoundary;

};

} //end namespace polytope

#endif
#endif
