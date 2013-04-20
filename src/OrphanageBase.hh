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


namespace polytope {

// Forward declaration
template<int Dimension, typename RealType> class Tessellator;


template<int Dimension, typename RealType>
class OrphanageBase
{
public:


  // Constructor
  OrphanageBase(const Tessellator<Dimension, RealType>* tessellator):
    mTessellator(tessellator) {
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
  void callPrivateTessellate(const std::vector<RealType>& points,
                             const std::vector<RealType>& PLCpoints,
                             const PLC<Dimension, RealType>& geometry,
                             const RealType* low,
                             const RealType* high,
                             const RealType dx,
                             Tessellation<Dimension, RealType>& mesh) const
  {
    mTessellator->tessellate(points, PLCpoints, geometry, low, high, dx, mesh);
  }

private:
  //-------------------- Private interface ----------------------//

  // Hold a pointer to a tessellator
  const Tessellator<Dimension, RealType>* mTessellator;
};

} // end polytope namespace

#endif
