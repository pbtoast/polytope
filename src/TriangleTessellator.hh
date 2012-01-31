//---------------------------------Spheral++----------------------------------//
// TriangleTessellator
// 
// An implemenation of the Tessellator interface that uses the Triangle
// library by Jonathan Shewchuk.
//----------------------------------------------------------------------------//
#ifndef __Polytope_TriangleTessellator__
#define __Polytope_TriangleTessellator__

#ifdef HAVE_TRIANGLE

#include <vector>
#include <cmath>

#include "Tessellator.hh"

namespace polytope 
{

template<typename Real>
class TriangleTessellator: public Tessellator<2, Real> 
{
  public:

  // Constructor, destructor.
  TriangleTessellator();
  ~TriangleTessellator();

  // Tessellate the given generators within the given bounding box.
  void tessellate(std::vector<Real>& points,
                  Real* low, Real* high,
                  Tessellation<2, Real>& mesh) const;

  // Tessellate the given generators. A bounding box is constructed about
  // the generators, and the corners of the bounding box are added as 
  // additional generators if they are not present in the list.
  void tessellate(std::vector<Real>& points,
                  Tessellation<2, Real>& mesh) const;

  // Tessellate obeying the given boundaries.
  void tessellate(std::vector<Real>& points,
                  const PLC<2, Real>& geometry,
                  Tessellation<2, Real>& mesh) const;

  // This Tessellator handles PLCs!
  bool handlesPLCs() const { return true; }

};

}

#endif
#endif
