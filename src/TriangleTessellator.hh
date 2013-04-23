
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
#include "Clipper2d.hh"
#include "BoostOrphanage.hh"
#include "Point.hh"
#include "polytope_tessellator_utilities.hh"
struct triangulateio;

// We use the Boost.Geometry library to handle polygon intersections and such.
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>

typedef int64_t CoordHash;
  
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

  // Return the tessellator name
  std::string name() const { return "TriangleTessellator"; }

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
  typedef boost::geometry::model::multi_polygon<BGpolygon> BGmulti_polygon;

  static CoordHash coordMax;
  static double degeneracy;

  // ------------------------------------------------- //
  // Specialized tessellations based on the point set  //
  // ------------------------------------------------- //

  // Compute node IDs around collinear generators and their quantized locations
  void computeCellNodesCollinear(const std::vector<RealType>& points,
				 std::map<IntPoint, std::pair<int,int> >& circumcenterMap,
				 std::map<int, std::vector<unsigned> >& cellNodes,
				 std::vector<unsigned>& infNodes) const;

  // Compute node IDs around each generator and their quantized locations
  void computeCellNodes(const std::vector<RealType>& points,
			std::map<IntPoint, std::pair<int,int> >& circumceterMap,
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

  // Bounding box used to quantize mesh nodes and mitigate degeneracies
  mutable std::vector<RealType> mLow, mHigh;
  mutable RealType mDelta;
   
  // Outer bounding box to quantize extreme circumcenters
  mutable std::vector<RealType> mLowOuter, mHighOuter;
  mutable RealType mDeltaOuter;

  // Infinite bounding circle
  mutable std::vector<RealType> mCenter;
  mutable RealType mRinf;

  // Maximum integer coordinate value
  mutable CoordHash mCoordMax;
  mutable RealType mDegeneracy;

  friend class BoostOrphanage<RealType>;

};

} //end polytope namespace

#endif
#endif
