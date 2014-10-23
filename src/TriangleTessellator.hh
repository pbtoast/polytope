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
#include "QuantTessellation.hh"

#include "Clipper2d.hh"
#include "BoostOrphanage.hh"
#include "QuantizedCoordinates.hh"
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

  // Tessellate obeying the given reducedPLC boundary.
  void tessellate(const std::vector<RealType>& points,
		  const ReducedPLC<2, RealType>& geometry,
		  Tessellation<2, RealType>& mesh) const;

  // This Tessellator handles PLCs!
  bool handlesPLCs() const { return true; }

  // Return the tessellator name
  std::string name() const { return "TriangleTessellator"; }

  //! Returns the accuracy to which this tessellator can distinguish coordinates.
  //! Should be returned appropriately for normalized coordinates, i.e., if all
  //! coordinates are in the range xi \in [0,1], what is the minimum allowed 
  //! delta in x.
  virtual RealType degeneracy() const { return mDegeneracy; }
  void degeneracy(RealType degeneracy) const { mDegeneracy = degeneracy; }

private:
  //-------------------- Private interface ---------------------- //
   
  // Typedefs for coordinates, edges and points
  typedef int64_t CoordHash;  
  typedef std::pair<int, int> EdgeHash;
  typedef Point2<CoordHash> IntPoint;
  typedef Point2<double> RealPoint;

  // Typedefs for boost.geometry rings and polygons
  // Template arguments are point_type and CW_orientation
  typedef boost::geometry::model::polygon<IntPoint ,false> BGpolygon;
  typedef boost::geometry::model::polygon<RealPoint,false> RealPolygon;
  typedef boost::geometry::model::ring<IntPoint, false>    BGring;
  typedef boost::geometry::model::ring<RealPoint,false>    RealRing;

  // ------------------------------------------------- //
  // Specialized tessellations based on the point set  //
  // ------------------------------------------------- //

  // Compute node IDs around each generator and their quantized locations
  void computeCellNodes(const std::vector<RealType>& points,
			std::vector<RealPoint>& nodeList,
			std::vector<std::vector<unsigned> >& cellNodes,
			std::vector<unsigned>& infNodes) const;
    
  // Compute bounded cell rings from collection of unbounded node locations
  void computeCellRings(const std::vector<RealType>& points,
			const std::vector<RealPoint>& nodeList,
			std::vector<std::vector<unsigned> >& cellNodes,
                        Clipper2d<CoordHash>& clipper,
			std::vector<BGring>& cellRings) const;

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
                  const std::vector<CoordHash>& IntPLCpoints,
                  const PLC<2, RealType>& geometry,
                  const QuantizedCoordinates<2, RealType>& coords,
                  std::vector<std::vector<std::vector<CoordHash> > >& IntCells) const;

  // -------------------------- //
  // Private member variables   //
  // -------------------------- //

  // The quantized coordinates for this tessellator (inner and outer)
  static CoordHash coordMax;
  static RealType mDegeneracy; 
  mutable QuantizedCoordinates<2,RealType> mCoords; //, mOuterCoords;

  friend class BoostOrphanage<RealType>;
};

} //end polytope namespace

#endif
#endif
