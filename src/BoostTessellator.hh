//------------------------------------------------------------------------
// BoostTessellator
// 
// Polytope wrapper for the native 2D Voronoi tessellator in Boost.Polygon
// v1.52 or greater
//------------------------------------------------------------------------
#ifndef __Polytope_BoostTessellator__
#define __Polytope_BoostTessellator__

#if HAVE_BOOST_VORONOI

#include <vector>
#include <cmath>

#include "Tessellator.hh"
#include "BoostOrphanage.hh"
#include "Point.hh"

// The Voronoi tools in Boost.Polygon
#include <boost/polygon/voronoi.hpp>

// We use the Boost.Geometry library to handle polygon intersections and such.
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>

// Some useful typedefs
typedef int64_t CoordHash;
typedef std::pair<int, int> EdgeHash;
typedef polytope::Point2<CoordHash> IntPoint;
typedef polytope::Point2<double> RealPoint;

typedef boost::polygon::voronoi_diagram<double> VD;

//------------------------------------------------------------------------
// Map Polytope's point class to Boost.Polygon
//------------------------------------------------------------------------
namespace boost{
namespace polygon{

template <>
struct geometry_concept<IntPoint> { typedef point_concept type; };
  
template <>
struct point_traits<IntPoint> {
  typedef CoordHash coordinate_type;
   
  static inline coordinate_type get(const IntPoint& point, orientation_2d orient) {
    return (orient == HORIZONTAL) ? point.x : point.y;
  }
};

} //end boost namespace
} //end polygon namespace


namespace polytope
{

template<typename RealType>
class BoostTessellator: public Tessellator<2, RealType>
{
public:

  // Constructor, destructor.
  BoostTessellator();
  ~BoostTessellator();

  // Tessellate the given generators. A bounding box is constructed about
  // the generators, and the corners of the bounding box are added as 
  // additional generators if they are not present in the list.
  void tessellate(const std::vector<RealType>& points,
                  Tessellation<2, RealType>& mesh) const;

  // Tessellate with a bounding box representing the boundaries.
  void tessellate(const std::vector<RealType>& points,
                  RealType* low,
                  RealType* high,
                  Tessellation<2, RealType>& mesh) const;

  // Tessellate obeying the given boundaries.
  void tessellate(const std::vector<RealType>& points,
                  const std::vector<RealType>& PLCpoints,
                  const PLC<2, RealType>& geometry,
                  Tessellation<2, RealType>& mesh) const;

  // This Tessellator handles PLCs!
  bool handlesPLCs() const { return true; }

  // The name of the tessellator
  std::string name() const { return "BoostTessellator"; }

private:
  //-------------------- Private interface ---------------------- //

  // Handy Boost.Geometry typedefs
  typedef boost::geometry::model::polygon<IntPoint,    // point type
                                          false>       // clockwise
    BGpolygon;
  typedef boost::geometry::model::ring<IntPoint,       // point type
                                       false>          // clockwise
    BGring;

  static CoordHash coordMax;
  static double degeneracy;

  // ------------------------------------------------- //
  // Specialized tessellations based on the point set  //
  // ------------------------------------------------- //

  // Compute nodes around collinear generators
  void computeCellNodesCollinear(const std::vector<RealType>& points,
                                 std::map<IntPoint, std::pair<int, int> >& nodeMap,
                                 std::map<int, std::vector<unsigned> >& cellNodes) const;

  // Compute the nodes around a generator
  void computeCellNodes(const std::vector<RealType>& points,
                         std::map<IntPoint, std::pair<int, int> >& nodeMap,
                         std::map<int, std::vector<unsigned> >& cellNodes) const;

  // Compute bounded cell rings from Boost Voronoi diagram
  void computeCellRings(const std::vector<RealType>& points,
			const std::vector<RealType>& PLCpoints,
			const PLC<2, RealType>& geometry,
			std::vector<BGring>& cellRings,
			bool performCellAdoption) const;

  // Compute an unbounded tessellation
  void computeVoronoiUnbounded(const std::vector<RealType>& points,
			       Tessellation<2, RealType>& mesh) const;

  // Compute a bounded tessellation
  void computeVoronoiBounded(const std::vector<RealType>& points,
			     const std::vector<RealType>& PLCpoints,
			     const PLC<2, RealType>& geometry,
			     Tessellation<2, RealType>& mesh) const;

  // ----------------------------------------------------- //
  // Private tessellate calls used by internal algorithms  //
  // ----------------------------------------------------- //

  // Bounded tessellation with prescribed bounding box
  void tessellate(const std::vector<RealType>& points,
                  const std::vector<RealType>& PLCpoints,
                  const PLC<2, RealType>& geometry,
                  const RealType* low,
                  const RealType* high,
                  const RealType dx,
                  Tessellation<2, RealType>& mesh) const;

  // -------------------------- //
  // Private member variables   //
  // -------------------------- //

  // Bounding box used to quantize points and mitigate degeneracies
  mutable std::vector<RealType> mLow, mHigh;
  mutable RealType mDelta;

  // Infinite bounding circle
  mutable std::vector<RealType> mCenter;
  mutable RealType mRinf;

  // Maximum integer coordinate value
  mutable CoordHash mCoordMax;

  friend class BoostOrphanage<RealType>;

};

} //end polytope namespace

#endif
#endif
