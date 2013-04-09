//------------------------------------------------------------------------
// TriangleTessellator
// 
// An implemenation of the Tessellator interface that uses the Triangle
// library by Jonathan Shewchuk.
//------------------------------------------------------------------------
#ifndef __Polytope_TriangleTessellator__
#define __Polytope_TriangleTessellator__

#if HAVE_TRIANGLE

#include <vector>
#include <cmath>

#include "Tessellator.hh"
#include "Point.hh"
struct triangulateio;

// We use the Boost.Geometry library to handle polygon intersections and such.
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/multi/geometries/multi_point.hpp>
#include <boost/geometry/multi/geometries/multi_polygon.hpp>


typedef int64_t CoordHash;
BOOST_GEOMETRY_REGISTER_POINT_2D(polytope::Point2<CoordHash>, CoordHash, boost::geometry::cs::cartesian, x, y);
BOOST_GEOMETRY_REGISTER_POINT_2D(polytope::Point2<double>, double, boost::geometry::cs::cartesian, x, y);

namespace polytope 
{

template<typename RealType>
class TriangleTessellator: public Tessellator<2, RealType> 
{
public:

  // Constructor, destructor.
  TriangleTessellator();
  ~TriangleTessellator();

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
   

private:
  //-------------------- Private interface ---------------------- //

  typedef std::pair<int, int> EdgeHash;
  typedef Point2<CoordHash> IntPoint;
  typedef Point2<double> RealPoint;
  typedef boost::geometry::model::polygon<IntPoint,    // point type
                                          false>       // clockwise
    BGpolygon;
  typedef boost::geometry::model::ring<IntPoint,       // point type
                                       false>          // clockwise
    BGring;
  typedef boost::geometry::model::polygon<RealPoint,   // point type
                                          false>       // clockwise
    realBGpolygon;
  typedef boost::geometry::model::multi_point<IntPoint> BGmulti_point;
  typedef boost::geometry::model::multi_polygon<BGpolygon> BGmulti_polygon;

  static CoordHash coordMax;
  static double degeneracy;

  // ------------------------------------------------- //
  // Specialized tessellations based on the point set  //
  // ------------------------------------------------- //

  // Compute an unbounded tessellation
  void computeVoronoiUnbounded(const std::vector<RealType>& points,
			       Tessellation<2, RealType>& mesh) const;

  // Compute a bounded tessellation
  void computeVoronoiBounded(const std::vector<RealType>& points,
			     const std::vector<RealType>& PLCpoints,
			     const PLC<2, RealType>& geometry,
			     Tessellation<2, RealType>& mesh) const;

  // Compute node IDs around each generator and their quantized locations
  void computeCellNodes(const std::vector<RealType>& points,
			std::map<IntPoint, std::pair<int,int> >& circumceterMap,
			std::map<int, std::vector<unsigned> >& cellNodes,
			std::vector<unsigned>& infNodes) const;
  
  
  // Compute node IDs around collinear generators and their quantized locations
  void computeCellNodesCollinear(const std::vector<RealType>& points,
				 std::map<IntPoint, std::pair<int,int> >& circumcenterMap,
				 std::map<int, std::vector<unsigned> >& cellNodes,
				 std::vector<unsigned>& infNodes) const;

  
  // Compute bounded cell rings from collection of unbounded node locations
  void computeCellRings(const std::vector<RealType>& points,
			const std::vector<RealType>& PLCpoints,
			const PLC<2, RealType>& geometry,
			std::vector<BGring>& cellRings,
			bool performCellAdoption) const;

  // Computes the triangularization using Triangle
  void computeDelaunay(const std::vector<RealType>& points,
                       triangulateio& delaunay) const;
  
  // Algorithm for assigning orphaned cells to neighboring cell rings
  void adoptionAlgorithm(const std::vector<RealType>& points,
			 std::vector<BGring>& cellRings,
			 std::vector<BGring>& orphanage) const;
  
  // Bounding box used to quantize points and mitigate degeneracies
  mutable std::vector<RealType> mLowOuter, mLowInner;
  mutable std::vector<RealType> mHighOuter, mHighInner;
  mutable RealType mdxOuter, mdxInner;

};

}

#endif
#endif
