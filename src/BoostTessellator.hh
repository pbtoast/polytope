//------------------------------------------------------------------------
// BoostTessellator
// 
// Polytope wrapper for the native 2D Voronoi tessellator in Boost.Polygon
// v1.52 or greater
//------------------------------------------------------------------------
#ifndef __Polytope_BoostTessellator__
#define __Polytope_BoostTessellator__

#ifdef HAVE_BOOST_VORONOI

#include <vector>
#include <cmath>
#include <limits>

#include "boost/polygon/voronoi.hpp"

#include "Tessellator.hh"
#include "QuantizedTessellation2d.hh"
#include "BoostOrphanage.hh"
#include "Point.hh"
#include "polytope_tessellator_utilities.hh"

namespace polytope {

template<typename RealType>
class BoostTessellator: public Tessellator<2, RealType> {
public:

  // The Boost.Polygon Voronoi diagram
  typedef boost::polygon::voronoi_diagram<RealType> VD;

  // Some useful typedefs
  typedef int                 CoordHash;
  typedef std::pair<int, int> EdgeHash;

  // The typedefs that follow from this choice
  typedef Point2<RealType>  RealPoint;
  typedef Point2<CoordHash> IntPoint;

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

  // Tessellate obeying the given PLC + points boundary.
  void tessellate(const std::vector<RealType>& points,
                  const std::vector<RealType>& PLCpoints,
                  const PLC<2, RealType>& geometry,
                  Tessellation<2, RealType>& mesh) const;

  // Tessellate obeying the given ReducedPLC boundary.
  void tessellate(const std::vector<RealType>& points,
                  const ReducedPLC<2, RealType>& geometry,
                  Tessellation<2, RealType>& mesh) const;

  // This Tessellator handles PLCs!
  bool handlesPLCs() const { return true; }

  // The name of the tessellator
  std::string name() const { return "BoostTessellator"; }

  //! Returns the accuracy to which this tessellator can distinguish coordinates.
  //! Should be returned appropriately for normalized coordinates, i.e., if all
  //! coordinates are in the range xi \in [0,1], what is the minimum allowed 
  //! delta in x.
  virtual RealType degeneracy() const { return mDegeneracy; }
  void degeneracy(const RealType val) const { mDegeneracy = val; }

private:
  //-------------------- Private interface ---------------------- //

  // ------------------------------------------------- //
  // Specialized tessellations based on the point set  //
  // ------------------------------------------------- //

  // Compute the nodes around a collection of generators
  void computeQuantTessellation(QuantizedTessellation2d<CoordHash, RealType>& result) const;

  // Compute the nodes around a linear, 1d collection of generators
  void computeCollinearQuantTessellation(QuantizedTessellation2d<CoordHash, RealType>& result) const;

  // void constructBoundedTopology(const std::vector<RealType>& points,
  //                               const ReducedPLC<2, RealType>& geometry,
  //                               const std::vector<ReducedPLC<2, CoordType> >& cellRings,
  //                               Tessellation<2, RealType>& mesh) const;

  // // ----------------------------------------------------- //
  // // Private tessellate calls used by internal algorithms  //
  // // ----------------------------------------------------- //

  // // Bounded tessellation with prescribed bounding box
  // void tessellate(const std::vector<RealType>& points,
  //                 const ReducedPLC<2, CoordHash>& intGeometry,
  //                 const QuantizedCoordinates<2,RealType>& coords,
  //                 std::vector<ReducedPLC<2, CoordHash> >& intCells) const;

  // -------------------------- //
  // Private member variables   //
  // -------------------------- //

  // The quantized coordinates for this tessellator (inner and outer)
  // static CoordHash coordMin, coordMax;
  // static RealType mxmin[2], mxmax[2], mlength;
  static RealType mDegeneracy; 

  friend class BoostOrphanage<RealType>;
};

//------------------------------------------------------------------------------
// Static initializations.
//------------------------------------------------------------------------------
// template<typename RealType> 
// typename BoostTessellator<RealType>::CoordHash
// BoostTessellator<RealType>::coordMin = std::numeric_limits<CoordHash>::min()/2;

// template<typename RealType> 
// typename BoostTessellator<RealType>::CoordHash
// BoostTessellator<RealType>::coordMax = std::numeric_limits<CoordHash>::max()/2;

template<typename RealType> 
RealType  
BoostTessellator<RealType>::mDegeneracy = 4.0/std::numeric_limits<CoordHash>::max();

} //end polytope namespace

#endif
#endif
