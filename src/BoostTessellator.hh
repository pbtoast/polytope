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
#include "Clipper2d.hh"
#include "BoostOrphanage.hh"
#include "Point.hh"
#include "polytope_tessellator_utilities.hh"

namespace polytope {

//------------------------------------------------------------------------------
// Define an intermediate struct to hold the quantized Voronoi.
//------------------------------------------------------------------------------
template<typename IntType, typename RealType>
struct QuantizedTessellation {
  typedef Point2<IntType> IntPoint;
  typedef Point2<RealType> RealPoint;
  
  RealType xmin[2], xmax[2], length, infRadius;
  std::vector<IntPoint> generators;

  std::vector<IntPoint> nodes;
  std::vector<std::pair<int, int> > edges;
  std::vector<RealPoint> infEdgeDirection;
  std::vector<std::vector<int> > cellEdges;

  // Construct with the given generators.  Finds the bounding limits, sets the infRadius,
  // and sets the quantized generators.
  QuantizedTessellation(const std::vector<RealType>& points) {
    geometry::computeBoundingBox<2, RealType>(&points[0], points.size(), true, xmin, xmax);
    length = std::max(result.xmax[0] - result.xmin[0], result.xmax[1] - result.xmin[1]);
    xmin[0] -= 2.0*length;
    xmin[1] -= 2.0*length;
    xmax[0] += 2.0*length;
    xmax[1] += 2.0*length;
    infRadius = 1.5*length;
    length *= 5.0;
    const int numGenerators = points.size()/2;
    generators.resize(numGenerators);
    for (unsigned i = 0; i < numGenerators; ++i) {
      this->quantize(&points[2*i], &generators[i].x);
      generators[i].index = i;
    }
  }

  // Convert real coordinates to integers.
  void quantize(const RealType* realcoords, IntType* intcoords) const {
    const RealType dx = length/(std::numeric_limits<IntType>::max() - std::numeric_limits<IntType>::min() + 1);
    intcoords[0] = std::numeric_limits<IntType>::min() + IntType((realcoords[0] - xmin[0])/dx);
    intcoords[1] = std::numeric_limits<IntType>::min() + IntType((realcoords[1] - xmin[1])/dx);
  }

  // Convert int coordinates to reals.
  void dequantize(const IntType* intcoords, RealType* realcoords) {
    const RealType dx = length/(std::numeric_limits<IntType>::max() - std::numeric_limits<IntType>::min() + 1);
    realcoords[0] = xmin[0] + (intcoords[0] - std::numeric_limits<IntType>::min())*dx;
    realcoords[1] = xmin[1] + (intcoords[1] - std::numeric_limits<IntType>::min())*dx;
  }
};

// The actual class.
  
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
  void computeQuantTessellation(const std::vector<RealType>& points,
                                QuantizedTessellation<CoordHash, RealType>& result) const;

  // Compute the nodes around a linear, 1d collection of generators
  void computeCellNodesCollinear(const std::vector<RealType>& points,
                                 QuantizedTessellation<CoordHash, RealType>& result) const;

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
// BoostTessellator<RealType>::coordMin = std::numeric_limits<CoordHash>::min();

// template<typename RealType> 
// typename BoostTessellator<RealType>::CoordHash
// BoostTessellator<RealType>::coordMax = std::numeric_limits<CoordHash>::max();

template<typename RealType> 
RealType  
BoostTessellator<RealType>::mDegeneracy = 1.0/std::numeric_limits<CoordHash>::max();

} //end polytope namespace

#endif
#endif
