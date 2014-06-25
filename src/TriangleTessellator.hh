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
#include "polytope_tessellator_utilities.hh"

struct triangulateio;

// We use the Boost.Geometry library to handle polygon intersections and such.
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>

namespace polytope 
{

template<typename RealType>
class TriangleTessellator: public Tessellator<2, RealType> 
{
public:

  typedef int64_t CoordHash;
  
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
  virtual RealType degeneracy() const { return 1.0e-8; }

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

  // ----------------------------------------------------- //
  // Private tessellate calls used by internal algorithms  //
  // ----------------------------------------------------- //
  // Computes the triangularization using Triangle
  void computeDelaunay(const std::vector<RealType>& points,
                       triangulateio& delaunay) const;

  // Internal method to compute the quantized tessellation.
  void
  computeUnboundedQuantizedTessellation(const std::vector<RealType>& points,
					const std::vector<RealType>& nonGeneratingPoints,
					internal::QuantTessellation<2, RealType>& qmesh) const;


  // Internal method to compute the quantized tessellation.
  void
  computeUnboundedQuantizedTessellation(const std::vector<RealType>& points,
					const std::vector<RealType>& nonGeneratingPoints,
					internal::QuantTessellation<2, RealType>& qmesh) const;

  
  // Internal method to compute the quantized tessellation.
  void
  computeDelaunayConnectivity(const std::vector<RealType>& points,
			      std::vector<RealPoint>& circumcenters,
			      std::vector<unsigned>& triMask,
			      std::map<EdgeHash, std::vector<unsigned> >& edge2tris,
			      std::map<int, std::set<unsigned> >& gen2tri,
			      std::vector<int>& triangleList,
			      RealPoint& low_inner,
			      RealPoint& high_inner,
			      RealPoint& low_outer,
			      RealPoint& high_outer) const;


  // -------------------------- //
  // Private member variables   //
  // -------------------------- //
  static CoordHash coordMax;
  static RealType mDegeneracy;

  friend class BoostOrphanage<RealType>;
};


template<typename RealType>
int64_t TriangleTessellator<RealType>::coordMax = (1LL << 26);

template<typename RealType>
RealType TriangleTessellator<RealType>::mDegeneracy = 1.0/TriangleTessellator::coordMax;



} //end polytope namespace

#endif
#endif
