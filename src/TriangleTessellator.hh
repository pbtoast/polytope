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

};

}

#endif
#endif
