//------------------------------------------------------------------------
// OrphanageBase
// 
// When intersecting unbounded cells with a PLC boundary, the resulting
// cell geometry may be in more than one piece. We keep the piece that
// contains the generator to be the new bounded cell. Any additional
// piece(s), which we term "orphans", are not included in the tessellation
//
// This is a base class for implementing algorithms that resolve the
// orphaned cell problem.
//
// Base class is templated on the dimension, real type, and a ring type.
// A ring is simply a vector of point type objects. For our purposes,
// ring points represent the vertices of a closed polygon.
//------------------------------------------------------------------------
#ifndef __Polytope_OrphanageBase__
#define __Polytope_OrphanageBase__

#include <set>
#include <map>
#include <vector>
#include "DimensionTraits.hh"
#include "QuantizedCoordinates.hh"


namespace polytope {

// Forward declaration
template<int Dimension, typename RealType> class Tessellator;


template<int Dimension, typename RealType>
class OrphanageBase
{
public:


  // Constructor
  OrphanageBase(const Tessellator<Dimension, RealType>* tessellatorPtr):
    mTessellatorPtr(tessellatorPtr) {
  };

  // Destructor
  virtual ~OrphanageBase() {};

  virtual void adoptOrphans(const std::vector<RealType>& points,
                            const RealType* low,
                            const RealType* high,
                            const RealType dx) const
  {
    error("This Orphanage does not support cell adoption");
  }

protected:

  // Call the underlying tessellator
  typedef typename DimensionTraits<Dimension, RealType>::CoordHash CoordHash;
  void callPrivateTessellate(const std::vector<RealType>& points,
                             const ReducedPLC<Dimension, CoordHash>& intGeometry,
                             const QuantizedCoordinates<Dimension, RealType>& coords,
                             std::vector<ReducedPLC<Dimension, CoordHash> >& intCells) const
  {
    mTessellatorPtr->tessellate(points, intGeometry, coords, intCells);
  }

private:
  //-------------------- Private interface ----------------------//

  // Hold a pointer to a tessellator
  const Tessellator<Dimension, RealType>* mTessellatorPtr;
};

} // end polytope namespace

#endif
